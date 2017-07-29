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
    nsearch merge --forward=<forward.fastq> --reverse=<reverse.fastq> --out=<merged.fastq>
    nsearch search --query=<query.fasta> --database=<database.fasta> --alnout=<output.aln> --minidentity=<minidentity> [--maxaccepts=<maxaccepts>] [--maxrejects=<maxrejects>]

  Options:
    --minidentity=<minidentity>    Minimum identity threshold (e.g. 0.8).
    --maxaccepts=<maxaccepts>      Maximum number of successful hits reported for one query [default: 1].
    --maxrejects=<maxrejects>      Abort after this many candidates were rejected [default: 8].
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

    Search( args[ "--query" ].asString(),
        args[ "--database" ].asString(),
        args[ "--alnout" ].asString(),
        std::stof( args[ "--minidentity" ].asString() ),
        args[ "--maxaccepts" ].asLong(),
        args[ "--maxrejects" ].asLong()
        );

    gStats.StopTimer();

    PrintSummaryHeader();
    PrintSummaryLine( gStats.ElapsedMillis() / 1000.0, "Seconds" );
  }

  // Merge
  if( args[ "merge" ].asBool() ) {
    gStats.StartTimer();

    Merge( args[ "--forward" ].asString(),
        args[ "--reverse" ].asString(),
        args[ "--out" ].asString() );

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
