#include <docopt.h>
#include <iostream>
#include <utility>
#include <functional>

#include <nsearch/Sequence.h>
#include <nsearch/FASTQ/Writer.h>
#include <nsearch/FASTA/Reader.h>
#include <nsearch/PairedEnd/Merger.h>
#include <nsearch/PairedEnd/Reader.h>
#include <nsearch/Database.h>

#include <nsearch/Alignment/ExtendAlign.h>
#include <nsearch/Alignment/BandedAlign.h>

#include "Stats.h"
#include "WorkerQueue.h"

Stats gStats;

static const char USAGE[] =
R"(
  Process and search sequences.

  Usage:
    nsearch merge <forward.fastq> <reverse.fastq> <merged.fastq>
    nsearch search <query.fasta> <database.fasta>
)";

void PrintProgressLine( const std::string label, size_t numProcessedReads, size_t numTotalReads )
{
  static std::string lastLabel;
  static int PROGRESS_BAR_WIDTH = 50;

  double done = double( numProcessedReads ) / double( numTotalReads );
  int pos = done * PROGRESS_BAR_WIDTH;

  if( label != lastLabel ) {
    std::cout << std::endl;
    lastLabel = label;
  }

  std::cout << "\r" << label << ": ";
  printf( "%.1f%%", done * 100.0 );
  std::cout << std::flush;
}

void PrintSummaryLine( double value, const char *line, double total = 0.0 )
{
  printf( "%15.1f  %s", value, line );
  if( total > 0.0 ) {
    printf( " (%.1f%%)", value / total * 100.0 );
  }
  printf( "\n" );
}

class QueuedWriter : public WorkerQueue< SequenceList > {
public:
  QueuedWriter( const std::string &path )
    : WorkerQueue( 1 ), mWriter( path )
  {
  }

protected:
  void Process( const SequenceList &list ) {
    for( auto seq : list ) {
      mWriter << seq;
    }
  }

private:
  FASTQ::Writer mWriter;
};

using PairedReads = std::pair< SequenceList, SequenceList >;

class QueuedMerger : public WorkerQueue< PairedReads > {
public:
  using ProcessedCallback = std::function< void ( size_t, size_t )  >;

  QueuedMerger( QueuedWriter &writer )
    : WorkerQueue( -1 ), mWriter( writer ), mTotalEnqueued( 0 ), mTotalProcessed( 0 )
  {
    mProcessedCallback = []( size_t, size_t ) { };
  }

  void Enqueue( PairedReads &item ) {
    mTotalEnqueued += item.first.size();
    WorkerQueue::Enqueue( item );
  }

  void OnProcessed( const ProcessedCallback &callback ) {
    mProcessedCallback = callback;
  }

protected:
  void Process( const std::pair< SequenceList, SequenceList > &queueItem ) {
    const SequenceList &fwd = queueItem.first;
    const SequenceList &rev = queueItem.second;

    const PairedEnd::Merger &merger = mMerger;

    Sequence mergedRead;
    SequenceList mergedReads;

    auto fit = fwd.begin();
    auto rit = rev.begin();
    while( fit != fwd.end() && rit != rev.end() ) {
      if( merger.Merge( mergedRead, *fit, *rit ) ) {
        gStats.numMerged++;
        gStats.mergedReadsTotalLength += mergedRead.Length();
        mergedReads.push_back( std::move( mergedRead ) );
      }

      mTotalProcessed++;
      gStats.numProcessed++;

      ++fit;
      ++rit;
    }

    if( !mergedReads.empty() )
      mWriter.Enqueue( mergedReads );

    mProcessedCallback( mTotalProcessed, mTotalEnqueued );
  }

private:
  size_t mTotalEnqueued;
  size_t mTotalProcessed;

  ProcessedCallback mProcessedCallback;
  QueuedWriter &mWriter;
  PairedEnd::Merger mMerger;
};

bool Merge( const std::string &fwdPath, const std::string &revPath, const std::string &mergedPath ) {
  PairedEnd::Reader reader( fwdPath, revPath );

  QueuedWriter writer( mergedPath  );
  QueuedMerger merger( writer );

  SequenceList fwdReads, revReads;

  while( !reader.EndOfFile() ) {
    reader.Read( fwdReads, revReads, 512 );
    auto pair = std::pair< SequenceList, SequenceList >( std::move( fwdReads ), std::move( revReads ) );
    merger.Enqueue( pair );
    PrintProgressLine( "Read Reads", reader.NumBytesRead(), reader.NumBytesTotal() );
  }

  merger.OnProcessed( []( size_t totalReadsProcessed, size_t totalReadsEnqueued ) {
    PrintProgressLine( "Merge Reads", totalReadsProcessed, totalReadsEnqueued );
  });

  merger.WaitTillDone();
  writer.WaitTillDone();

  return true;
}

bool Search( const std::string &queryPath, const std::string &databasePath ) {
  Sequence seq;
  SequenceList sequences;

  FASTA::Reader dbReader( databasePath );
  std::cout << "Indexing DB" << std::endl;
  while( !dbReader.EndOfFile() ) {
    dbReader >> seq;
    sequences.push_back( std::move( seq ) );
  }

  Database db( sequences, 8 );
  /* db.Stats(); */
  /* return false; */

  std::cout << "Querying DB" << std::endl;
  FASTA::Reader qryReader( queryPath );
  while( !qryReader.EndOfFile() )  {
    qryReader >> seq;
    SequenceList candidates = db.Query( seq, 0.75, 1, 8 );
    /* std::cout << seq.identifier << std::endl; */
    /* for( auto &candidate : candidates ) { */
    /*   std::cout << " " << candidate.identifier << std::endl; */
    /* } */
    /* std::cout << "===" << std::endl; */
  }

  return true;
}

int main( int argc, const char **argv ) {
  std::map<std::string, docopt::value> args
    = docopt::docopt(USAGE,
        { argv + 1, argv + argc },
        true, // help
        "nsearch");

  if( args[ "search" ].asBool() ) {
    gStats.StartTimer();

    Search( args[ "<query.fasta>" ].asString(),
        args[ "<database.fasta>" ].asString() );

    gStats.StopTimer();

    std::cout << std::endl;
    std::cout << "Summary:" << std::endl;
    PrintSummaryLine( gStats.ElapsedMillis() / 1000.0, "Seconds" );
  }

  if( args[ "merge" ].asBool() ) {
    gStats.StartTimer();

    Merge( args[ "<forward.fastq>" ].asString(),
        args[ "<reverse.fastq>" ].asString(),
        args[ "<merged.fastq>" ].asString() );


    gStats.StopTimer();

    std::cout << std::endl;
    std::cout << "Summary:" << std::endl;
    PrintSummaryLine( gStats.ElapsedMillis() / 1000.0, "Seconds" );
    PrintSummaryLine( gStats.numProcessed / gStats.ElapsedMillis(), "Processed/ms" );
    PrintSummaryLine( gStats.numProcessed, "Pairs" );
    PrintSummaryLine( gStats.numMerged, "Merged", gStats.numProcessed );
    PrintSummaryLine( gStats.MeanMergedLength(), "Mean merged length" );
  }

  return 0;
}
