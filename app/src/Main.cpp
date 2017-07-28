#include <docopt.h>
#include <iostream>
#include <utility>
#include <functional>
#include <sstream>

#include <nsearch/Sequence.h>
#include <nsearch/FASTA/Reader.h>

#include "Common.h"
#include "Stats.h"
#include "Merge.h"
#include "Search.h"

Stats gStats;

static const char USAGE[] =
R"(
  Process and search sequences.

  Usage:
    nsearch merge <forward.fastq> <reverse.fastq> <merged.fastq>
    nsearch search <query.fasta> <database.fasta> <output.txt>
)";


void PrintSummaryHeader()
{
  std::cout << std::endl << std::endl;
  std::cout << "Summary:" << std::endl;
}

void PrintSummaryLine( float value, const std::string &line, float total = 0.0, UnitType unit = UnitType::COUNTS )
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

int main( int argc, const char **argv ) {
  std::map<std::string, docopt::value> args
    = docopt::docopt( USAGE,
        { argv + 1, argv + argc },
        true, // help
        APP_NAME );

  // Print header
  std::cout << APP_NAME << " " << APP_VERSION << " (built on " << BUILD_TIMESTAMP << ")" << std::endl;

  // Search
  if( args[ "search" ].asBool() ) {
    gStats.StartTimer();

    Search( args[ "<query.fasta>" ].asString(),
        args[ "<database.fasta>" ].asString(),
        args[ "<output.txt>" ].asString() );

    gStats.StopTimer();

    PrintSummaryHeader();
    PrintSummaryLine( gStats.ElapsedMillis() / 1000.0, "Seconds" );
  }

  // Merge
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
