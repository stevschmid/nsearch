#include <docopt.h>
#include <functional>
#include <iostream>
#include <sstream>
#include <utility>

#include <nsearch/FASTA/Reader.h>
#include <nsearch/Sequence.h>
#include <nsearch/Alphabet/DNA.h>
#include <nsearch/Alphabet/Protein.h>

#include "Common.h"
#include "Filter.h"
#include "Merge.h"
#include "Search.h"
#include "Stats.h"

Stats gStats;

static const char USAGE[] = R"(
  Metagenomics tool for the rest of us.

  Usage:
    nsearch search --query=<queryfile> --db=<databasefile>
      --out=<outputfile> --min-identity=<minidentity> [--max-hits=<maxaccepts>] [--max-rejects=<maxrejects>] [--protein] [--strand=<strand>]
    nsearch merge --forward=<forwardfile> --reverse=<reversefile> --out=<outputfile>
    nsearch filter --in=<inputfile> --out=<outputfile> [--max-expected-errors=<maxee>]

  Options:
    --min-identity=<minidentity>    Minimum identity threshold (e.g. 0.8).
    --max-hits=<maxaccepts>         Maximum number of successful hits reported for one query [default: 1].
    --max-rejects=<maxrejects>      Abort after this many candidates were rejected [default: 16].
    --max-expected-errors=<maxee>   Maximum number of expected errors [default: 1.0].
    --strand=<strand>               Strand to search on (plus, minus or both). If minus (or both), queries are reverse complemented [default: both].
)";

void PrintSummaryHeader() {
  std::cout << std::endl << std::endl;
  std::cout << "Summary:" << std::endl;
}

void PrintSummaryLine( const float value, const std::string& line,
                       const float    total = 0.0,
                       const UnitType unit  = UnitType::COUNTS ) {
  std::ios::fmtflags f( std::cout.flags() );
  std::cout << std::setw( 10 ) << ValueWithUnit( value, unit );
  std::cout << ' ' << line;
  if( total > 0.0 ) {
    std::cout << " (" << value / total * 100.0 << "%)";
  }
  std::cout << std::endl;
  std::cout.flags( f );
}

using Args = std::map< std::string, docopt::value > ;

template < typename A >
void AddSpecialSearchParams( const Args& args, SearchParams< A >* sp ) {}

void AddSpecialSearchParams( const Args& args, SearchParams< DNA >* sp ) {
  auto str = args.at( "--strand" ).asString();
  sp->strand = DNA::Strand::Plus;
  if( str == "minus" ) {
    sp->strand = DNA::Strand::Minus;
  } else if( str == "both" ) {
    sp->strand = DNA::Strand::Both;
  }
}

template < typename A >
SearchParams< A > ParseSearchParams( const Args& args ) {
  SearchParams< A > sp;

  sp.minIdentity = std::stof( args.at( "--min-identity" ).asString() );
  sp.maxAccepts  = args.at( "--max-hits" ).asLong();
  sp.maxRejects  = args.at( "--max-rejects" ).asLong();

  AddSpecialSearchParams( args, &sp );

  return sp;
}

int main( int argc, const char** argv ) {
  Args args = docopt::docopt( USAGE, { argv + 1, argv + argc },
                              true, // help
                              APP_NAME );

  // Print header
  std::cout << APP_NAME << " " << APP_VERSION << " (built on "
            << BUILD_TIMESTAMP << ")" << std::endl;

  // Search
  if( args[ "search" ].asBool() ) {
    gStats.StartTimer();

    auto query      = args[ "--query" ].asString();
    auto db         = args[ "--db" ].asString();
    auto out        = args[ "--out" ].asString();


    if( args[ "--protein" ].asBool() ) {
      DoSearch< Protein >( query, db, out, ParseSearchParams< Protein >( args ) );
    } else {
      DoSearch< DNA >( query, db, out, ParseSearchParams< DNA >( args ) );
    }

    gStats.StopTimer();

    PrintSummaryHeader();
    PrintSummaryLine( gStats.ElapsedMillis() / 1000.0, "Seconds" );
  }

  // Merge
  if( args[ "merge" ].asBool() ) {
    gStats.StartTimer();

    DoMerge( args[ "--forward" ].asString(), args[ "--reverse" ].asString(),
           args[ "--out" ].asString() );

    gStats.StopTimer();

    PrintSummaryHeader();
    PrintSummaryLine( gStats.ElapsedMillis() / 1000.0, "Seconds" );
    PrintSummaryLine( gStats.numProcessed / gStats.ElapsedMillis(),
                      "Processed/ms" );
    PrintSummaryLine( gStats.numProcessed, "Pairs" );
    PrintSummaryLine( gStats.numMerged, "Merged", gStats.numProcessed );
    PrintSummaryLine( gStats.MeanMergedLength(), "Mean merged length" );
  }

  // Filter
  if( args[ "filter" ].asBool() ) {
    gStats.StartTimer();

    auto in    = args[ "--in" ].asString();
    auto out   = args[ "--out" ].asString();
    auto maxee = std::stof( args[ "--max-expected-errors" ].asString() );

    DoFilter( in, out, maxee );

    gStats.StopTimer();

    PrintSummaryHeader();
    PrintSummaryLine( gStats.ElapsedMillis() / 1000.0, "Seconds" );
  }

  return 0;
}
