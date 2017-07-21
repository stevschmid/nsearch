#include <docopt.h>
#include <iostream>
#include <utility>
#include <functional>
#include <sstream>

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

enum class UnitType { COUNTS, BYTES };

std::string ValueWithUnit( double value, UnitType unit ) {
  static const std::map< UnitType, std::map< size_t, std::string > > CONVERSION = {
    {
      UnitType::COUNTS,
      {
        { 1, "" },
        { 1000, "k" },
        { 1000 * 1000, "M" },
        { 1000 * 1000 * 1000, "G" }
      },
    },
    {
      UnitType::BYTES,
      {
        { 1, " Bytes" },
        { 1024, " kB" },
        { 1024 * 1024, " MB" },
        { 1024 * 1024 * 1024, " GB" }
      }
    }
  };

  auto list = CONVERSION.find( unit  )->second;
  auto it = list.begin();
  while( it != list.end() && value > it->first * 10 ) {
    it++;
  }

  std::stringstream ss;
  if( it != list.begin() ) {
    it--;
    value = floor( value / it->first );
  }

  if( it->first == 1 ) {
    ss << std::setprecision(1) << std::setiosflags( std::ios::fixed );
  }

  ss << value;
  ss << it->second;
  return ss.str();
}

void PrintProgressLine( const std::string label, size_t value, size_t max, UnitType unit = UnitType::COUNTS )
{
  static std::string lastLabel;
  static int PROGRESS_BAR_WIDTH = 50;

  std::ios::fmtflags f( std::cout.flags() );

  double done = double( value ) / double( max );
  int pos = done * PROGRESS_BAR_WIDTH;

  if( label != lastLabel ) {
    std::cout << std::endl;
    lastLabel = label;
  }

  std::cout << label << ": ";
  std::cout << done * 100.0 << '%';
  std::cout << " (" << ValueWithUnit( value, unit ) << ")";
  std::cout << std::string( 20, ' ' ) << "\r" << std::flush;
  std::cout.flags( f );
}

void PrintSummaryHeader()
{
  std::cout << std::endl << std::endl;
  std::cout << "Summary:" << std::endl;
}

void PrintSummaryLine( double value, const std::string &line, double total = 0.0, UnitType unit = UnitType::COUNTS )
{
  std::ios::fmtflags f( std::cout.flags() );
  std::cout << std::setw( 10 ) << ValueWithUnit( value, unit );
  std::cout << ' ' << line;
  if( total > 0.0 ) {
    std::cout << " (" << value / total * 100.0 << "%)";
  }
  std::cout << std::endl;
  std::cout.flags( f );
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
  QueuedMerger( QueuedWriter &writer )
    : WorkerQueue( -1 ), mWriter( writer )
  {
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

      ++fit;
      ++rit;
    }

    if( !mergedReads.empty() )
      mWriter.Enqueue( mergedReads );

    gStats.numProcessed += fwd.size();
  }

private:
  QueuedWriter &mWriter;
  PairedEnd::Merger mMerger;
};

bool Merge( const std::string &fwdPath, const std::string &revPath, const std::string &mergedPath ) {
  const int numReadsPerWorkItem = 512;

  PairedEnd::Reader reader( fwdPath, revPath );

  QueuedWriter writer( mergedPath  );
  QueuedMerger merger( writer );

  SequenceList fwdReads, revReads;

  while( !reader.EndOfFile() ) {
    reader.Read( fwdReads, revReads, numReadsPerWorkItem );
    auto pair = std::pair< SequenceList, SequenceList >( std::move( fwdReads ), std::move( revReads ) );
    merger.Enqueue( pair );
    PrintProgressLine( "Reading File", reader.NumBytesRead(), reader.NumBytesTotal(), UnitType::BYTES );
  }

  merger.OnProcessed( []( size_t totalReadsProcessed, size_t totalReadsEnqueued ) {
    PrintProgressLine( "Merging Reads", totalReadsProcessed * numReadsPerWorkItem, totalReadsEnqueued * numReadsPerWorkItem );
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
    = docopt::docopt( USAGE,
        { argv + 1, argv + argc },
        true, // help
        APP_NAME );

  // Print header
  std::cout << APP_NAME << " " << APP_VERSION << " (built on " << BUILD_TIMESTAMP << ")" << std::endl;

  // Show one decimal point
  std::cout << std::setiosflags( std::ios::fixed ) << std::setprecision( 1 );

  if( args[ "search" ].asBool() ) {
    gStats.StartTimer();

    Search( args[ "<query.fasta>" ].asString(),
        args[ "<database.fasta>" ].asString() );

    gStats.StopTimer();

    PrintSummaryHeader();
    PrintSummaryLine( gStats.ElapsedMillis() / 1000.0, "Seconds" );
  }

  if( args[ "merge" ].asBool() ) {
    gStats.StartTimer();

    Merge( args[ "<forward.fastq>" ].asString(),
        args[ "<reverse.fastq>" ].asString(),
        args[ "<merged.fastq>" ].asString() );


    gStats.StopTimer();

    PrintSummaryHeader();
    PrintSummaryLine( gStats.ElapsedMillis() / 1000.0, "Seconds" );
    PrintSummaryLine( gStats.numProcessed / gStats.ElapsedMillis(), "Processed/ms" );
    PrintSummaryLine( gStats.numProcessed, "Pairs" );
    PrintSummaryLine( gStats.numMerged, "Merged", gStats.numProcessed );
    PrintSummaryLine( gStats.MeanMergedLength(), "Mean merged length" );
  }

  return 0;
}
