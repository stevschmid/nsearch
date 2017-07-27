#include <docopt.h>
#include <iostream>
#include <utility>
#include <functional>
#include <sstream>

#include <nsearch/Sequence.h>
#include <nsearch/FASTA/Reader.h>

#include <nsearch/Database.h>
#include <nsearch/Database/GlobalSearch.h>

#include "Common.h"
#include "Stats.h"
#include "Merge.h"

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

bool Search( const std::string &queryPath, const std::string &databasePath, const std::string &outputPath ) {
  ProgressOutput progress;

  Sequence seq;
  SequenceList sequences;

  FASTA::Reader dbReader( databasePath );

  enum ProgressType { ReadDBFile, StatsDB, IndexDB, ReadQueryFile, SearchDB };

  progress.Add( ProgressType::ReadDBFile, "Reading DB file", UnitType::BYTES );
  progress.Add( ProgressType::StatsDB, "Analyzing DB sequences");
  progress.Add( ProgressType::IndexDB, "Indexing database");
  progress.Add( ProgressType::ReadQueryFile, "Reading query file", UnitType::BYTES );
  progress.Add( ProgressType::SearchDB, "Searching database" );

  // Read DB
  progress.Activate( ProgressType::ReadDBFile );
  while( !dbReader.EndOfFile() ) {
    dbReader >> seq;
    sequences.push_back( std::move( seq ) );
    progress.Set( ProgressType::ReadDBFile, dbReader.NumBytesRead(), dbReader.NumBytesTotal() );
  }

  // Index DB
  auto dbCallback = [&]( Database::ProgressType type, size_t num, size_t total ) {
    switch( type ) {
      case Database::ProgressType::StatsCollection:
        progress.Activate( ProgressType::StatsDB ).Set( ProgressType::StatsDB, num, total );
        break;

      case Database::ProgressType::Indexing:
        progress.Activate( ProgressType::IndexDB ).Set( ProgressType::IndexDB, num, total );
        break;

      default:
        break;
    }
  };
  Database db( sequences, 8, dbCallback );

  // Read query
  FASTA::Reader qryReader( queryPath );

  SequenceList queries;
  progress.Activate( ProgressType::ReadQueryFile );
  while( !qryReader.EndOfFile() )  {
    qryReader >> seq;
    progress.Set( ProgressType::ReadQueryFile, qryReader.NumBytesRead(), qryReader.NumBytesTotal() );
    queries.push_back( seq );
  }

  // Search
  progress.Activate( ProgressType::SearchDB );
  GlobalSearch search( db, 0.75, 1, 8 );

  size_t count = 0;

  std::ofstream of;
  of.open( outputPath );

  for( auto &query : queries ) {
    if( (count++) % 50 == 0 || count == queries.size() ) {
      progress.Set( ProgressType::SearchDB, count, queries.size() );
    }

    Search::ResultList results = search.Query( query );
    if( results.empty() )
      continue;

    of << "Found: " << results.size() << " hits" << std::endl << std::endl;
    for( auto &result : results ) {
      of << result << std::endl;
    }
    of << std::string( 50, '-') << std::endl << std::endl;
  }

  of.close();

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
