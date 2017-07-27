#include "Search.h"

#include <nsearch/Sequence.h>
#include <nsearch/FASTA/Reader.h>
#include <nsearch/Database.h>
#include <nsearch/Database/GlobalSearch.h>

#include "Common.h"

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

    auto results = search.Query( query );
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
