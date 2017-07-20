#include <docopt.h>
#include <iostream>
#include <utility>

#include <nsearch/Sequence.h>
#include <nsearch/FASTQ/Writer.h>
#include <nsearch/FASTA/Reader.h>
#include <nsearch/PairedEnd/Merger.h>
#include <nsearch/PairedEnd/Reader.h>
#include <nsearch/Database.h>
#include <nsearch/Aligner.h>

#include <nsearch/Alignment/ExtendAlign.h>
#include <nsearch/Alignment/BandedAlign.h>

#include "Stats.h"
#include "WorkerQueue.h"

#include <nsearch/SpacedSeeds.h>

Stats gStats;

static const char USAGE[] =
R"(
  Process and search sequences.

  Usage:
    nsearch merge <forward.fastq> <reverse.fastq> <merged.fastq>
    nsearch search <query.fasta> <database.fasta>
)";

void PrintProgressLine( size_t numProcessedReads, size_t numTotalReads )
{
  static int PROGRESS_BAR_WIDTH = 50;

  double done = double( numProcessedReads ) / double( numTotalReads );
  int pos = done * PROGRESS_BAR_WIDTH;

  std::cout << "\r[";
  std::cout << std::string( pos, '=' );
  if( done < 1.0 )
    std::cout << '>';
  std::cout << std::string( std::max( PROGRESS_BAR_WIDTH - pos - 1, 0 ), ' ' );
  std::cout << "] ";
  printf( " %.1f%%", done * 100.0 );
  std::cout << std::flush;
}

void PrintSummaryLine( double value, const char *line, double total = 0.0 )
{
  printf( "%10.1f  %s", value, line );
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

class QueuedMerger : public WorkerQueue< std::pair< SequenceList, SequenceList > > {
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

      gStats.numProcessed++;

      ++fit;
      ++rit;
    }

    if( !mergedReads.empty() )
      mWriter.Enqueue( mergedReads );
  }

private:
  QueuedWriter &mWriter;
  PairedEnd::Merger mMerger;
};

bool Merge( const std::string &fwdPath, const std::string &revPath, const std::string &mergedPath ) {
  PairedEnd::Reader reader( fwdPath, revPath );

  QueuedWriter writer( mergedPath );
  QueuedMerger merger( writer );

  SequenceList fwdReads, revReads;

  while( !reader.EndOfFile() ) {
    reader.Read( fwdReads, revReads, 512 );
    auto pair = std::pair< SequenceList, SequenceList >( std::move( fwdReads ), std::move( revReads ) );
    merger.Enqueue( pair );
    PrintProgressLine( reader.NumBytesRead(), reader.NumBytesTotal() );
  }

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
    SequenceList candidates = db.Query( seq, 0.75, 1, 20 );
    /* std::cout << seq.identifier << std::endl; */
    /* for( auto &candidate : candidates ) { */
    /*   std::cout << " " << candidate.identifier << std::endl; */
    /* } */
    /* std::cout << "===" << std::endl; */
  }

  return true;
}

/* template<typename T> */
/* void printBin(const T& a) */
/* { */
/*   const char* beg = reinterpret_cast<const char*>(&a); */
/*   const char* end = beg + sizeof(a); */
/*   while(beg != end) */
/*     std::cout << std::bitset<CHAR_BIT>(*beg++) << ' '; */
/*   std::cout << '\n'; */
/* } */

int main( int argc, const char **argv ) {
  std::map<std::string, docopt::value> args
    = docopt::docopt(USAGE,
        { argv + 1, argv + argc },
        true, // help
        "nsearch");

  // Test cases
  // A>>>>B
  // B>>>>A
  // Empty A
  // empty B
  // A breaking case when first row is not initialized properly (beyond bandwidth)
  // THIS CASE:
    /* Sequence A = "AAAAAAAAAAAAAAA"; */
    /* Sequence B = "CCCCCCAAAAAAAAA"; */
    /* int score = ba.Align( A, B, &cig, 0, 0, AlignmentDirection::forwards ); */
    /* std::cout << score << std::endl; */
    /* std::cout << cig << std::endl; */
    /* std::cout << A.sequence << std::endl; */
    /* A = "CCCCCCCCCCCCCCC"; */
    /* B = "CCCCCCAAAAAAAAA"; */
    /* score = ba.Align( A, B, &cig, 0, 0, AlignmentDirection::forwards ); */
    /* std::cout << score << std::endl; */
    /* std::cout << cig << std::endl; */
    /* std::cout << A.sequence << std::endl; */

  // TEST CASE FOR WHEN STARTA>>>LENA
    /* BandedAlignParams bap; */
    /* BandedAlign ba( bap ); */
    /* Cigar cig; */
    /* Sequence A = "ATGCC"; */
    /* Sequence B = "XXXATGCC"; */
    /* int score = ba.Align( A, B, &cig, AlignmentDirection::forwards, 6, 3 ); */
    /* std::cout << score << std::endl; */
    /* std::cout << cig << std::endl; */
    /* std::cout << A.sequence << std::endl; */
  /* return 0; */

  /* SpacedSeeds seedsC( "TUTGT", 4 ); */
  /* size_t c1 = 0; */
  /* seedsC.ForEach( [&]( size_t pos, uint32_t word ) { */
  /*   printBin( word ); */
  /*   std::cout << std::endl; */
  /*   c1++; */
  /* }); */

  /* std::cout << c1 << std::endl; */
  /* std::cout << "====" << std::endl; */

  /* SpacedSeedsSIMD seeds( "TUTGT", 4 ); */
  /* size_t c2 = 0; */
  /* seeds.ForEach( [&]( size_t pos, uint32_t word ) { */
  /*   printBin( word ); */
  /*   std::cout << std::endl; */
  /*   c2++; */
  /* }); */
  /* std::cout << c2 << std::endl; */
  /* return 0; */

  /* ExtendAlign ea; */
  /* ExtendedAlignment aln; */
  /* Sequence A = "AATTT"; */
  /* Sequence B = "GGGGT"; */
  /* size_t bestA, bestB; */
  /* Cigar cigar; */
  /* int score = ea.Extend( A, B, &bestA, &bestB, &cigar, AlignmentDirection::backwards, A.Length(), B.Length() ); */
  /* std::cout << score << std::endl; */
  /* std::cout << aln.cigar << std::endl; */
  /* std::cout << bestA << ", " << bestB << std::endl; */

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
