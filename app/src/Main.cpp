#include <docopt.h>
#include <iostream>

#include <nsearch/Sequence.h>
#include <nsearch/FASTQ/Writer.h>
#include <nsearch/PairedEnd/Merger.h>
#include <nsearch/PairedEnd/Reader.h>

#include "Stats.h"
#include "ThreadPool.h"

Stats gStats;

static const char USAGE[] =
R"(
  Process and search sequences.

  Usage:
    nsearch merge <forward.fastq> <reverse.fastq> <merged.fastq>
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

bool Merge( const std::string &fwdPath, const std::string &revPath, const std::string &mergedPath ) {
  PairedEnd::Merger mergerObj;
  const PairedEnd::Merger &merger = mergerObj; // const so we ensure thread-safety

  PairedEnd::Reader reader( fwdPath, revPath );
  FASTQ::Writer writer( mergedPath );

  ThreadPool pool;
  SequenceList fwdReads, revReads;

  while( !reader.EndOfFile() ) {
    reader.Read( fwdReads, revReads, 512 );

    auto mergeAndWrite = [ &merger, &writer ] ( SequenceList &fwd, SequenceList &rev ) {
      Sequence mergedRead;

      auto fit = fwd.begin();
      auto rit = rev.begin();
      while( fit != fwd.end() && rit != rev.end() ) {
        if( merger.Merge( mergedRead, *fit, *rit ) ) {
          gStats.numMerged++;
          gStats.mergedReadsTotalLength += mergedRead.Length();

          /* writer << mergedRead; */
        }

        gStats.numProcessed++;

        ++fit;
        ++rit;
      }
    };

    auto task = std::bind( mergeAndWrite, std::move( fwdReads ), std::move( revReads ) );
    pool.Enqueue( task );
    PrintProgressLine( reader.NumBytesRead(), reader.NumBytesTotal() );
  }

  while( !pool.Done() )
    std::this_thread::sleep_for( std::chrono::milliseconds( 50 ) );

  return true;
}

int main( int argc, const char **argv ) {
  std::map<std::string, docopt::value> args
    = docopt::docopt(USAGE,
        { argv + 1, argv + argc },
        true, // help
        "nsearch");

  if( args["merge"] ) {
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
