#include "Search.h"

#include <nsearch/Alnout/Writer.h>
#include <nsearch/Database.h>
#include <nsearch/Database/GlobalSearch.h>
#include <nsearch/FASTA/Reader.h>
#include <nsearch/Sequence.h>
#include <nsearch/Sequence/Protein.h>

#include "WorkerQueue.h"

#include "Common.h"

template < typename A >
using QueryWithHits     = std::pair< Sequence< A >, HitList< A > >;

template < typename A >
using QueryWithHitsList = std::deque< QueryWithHits< A > >;

template < typename A >
class QueueItemInfo< QueryWithHitsList< A > > {
public:
  static size_t Count( const QueryWithHitsList< A >& list ) {
    return std::accumulate(
      list.begin(), list.end(), 0,
      []( int sum, const QueryWithHits< A >& q ) { return sum + q.second.size(); } );
  }
};

template < typename A >
class SearchResultsWriterWorker {
public:
  SearchResultsWriterWorker( const std::string& path ) : mWriter( path ) {}

  void Process( const QueryWithHitsList< A >& queryWithHitsList ) {
    for( auto& queryWithHits : queryWithHitsList ) {
      mWriter << queryWithHits;
    }
  }

private:
  Alnout::Writer< A > mWriter;
};

template < typename A >
using SearchResultsWriter =
  WorkerQueue< SearchResultsWriterWorker< A >, QueryWithHitsList< A >,
               const std::string& >;

template < typename A >
class QueueItemInfo< SequenceList< A > > {
public:
  static size_t Count( const SequenceList< A >& list ) {
    return list.size();
  }
};

template < typename A >
class QueryDatabaseSearcherWorker {
public:
  QueryDatabaseSearcherWorker( SearchResultsWriter< A >* writer,
                               const Database< A >*      database,
                               const float minIdentity, const int maxAccepts,
                               const int maxRejects )
      : mWriter( *writer ),
        mGlobalSearch( *database, minIdentity, maxAccepts, maxRejects ) {}

  void Process( const SequenceList< A >& queries ) {
    QueryWithHitsList< A > list;

    for( auto& query : queries ) {
      auto hits = mGlobalSearch.Query( query );
      if( hits.empty() )
        continue;

      list.push_back( { query, hits } );
    }

    if( !list.empty() ) {
      mWriter.Enqueue( list );
    }
  }

private:
  GlobalSearch< A >         mGlobalSearch;
  SearchResultsWriter< A >& mWriter;
};

template < typename A >
using QueryDatabaseSearcher =
  WorkerQueue< QueryDatabaseSearcherWorker< A >, SequenceList< A >,
               SearchResultsWriter< A >*, const Database< A >*, const float,
               const int, const int >;

bool Search( const std::string& queryPath, const std::string& databasePath,
             const std::string& outputPath, const float minIdentity,
             const int maxAccepts, const int maxRejects ) {
  ProgressOutput progress;

  Sequence< Protein >     seq;
  SequenceList< Protein > sequences;

  FASTA::Reader< Protein > dbReader( databasePath );

  enum ProgressType {
    ReadDBFile,
    StatsDB,
    IndexDB,
    ReadQueryFile,
    SearchDB,
    WriteHits
  };

  progress.Add( ProgressType::ReadDBFile, "Read database", UnitType::BYTES );
  progress.Add( ProgressType::StatsDB, "Analyze database" );
  progress.Add( ProgressType::IndexDB, "Index database" );
  progress.Add( ProgressType::ReadQueryFile, "Read queries", UnitType::BYTES);
  progress.Add( ProgressType::SearchDB, "Search database" );
  progress.Add( ProgressType::WriteHits, "Write hits" );

  // Read DB
  progress.Activate( ProgressType::ReadDBFile );
  while( !dbReader.EndOfFile() ) {
    dbReader >> seq;
    sequences.push_back( std::move( seq ) );
    progress.Set( ProgressType::ReadDBFile, dbReader.NumBytesRead(),
                  dbReader.NumBytesTotal() );
  }

  // Index DB
  const int wordSize = 5;
  Database< Protein > db( wordSize );
  db.SetProgressCallback( [&]( Database< Protein >::ProgressType type, size_t num,
                         size_t total ) {
    switch( type ) {
      case Database< Protein >::ProgressType::StatsCollection:
        progress.Activate( ProgressType::StatsDB )
          .Set( ProgressType::StatsDB, num, total );
        break;

      case Database< Protein >::ProgressType::Indexing:
        progress.Activate( ProgressType::IndexDB )
          .Set( ProgressType::IndexDB, num, total );
        break;

      default:
        break;
    }
  });
  db.Initialize( sequences );

  // Read and process queries
  const int numQueriesPerWorkItem = 64;

  SearchResultsWriter< Protein >   writer( 1, outputPath );
  QueryDatabaseSearcher< Protein > searcher( -1, &writer, &db, minIdentity,
                                             maxAccepts, maxRejects );

  searcher.OnProcessed( [&]( size_t numProcessed, size_t numEnqueued ) {
    progress.Set( ProgressType::SearchDB, numProcessed, numEnqueued );
  } );
  writer.OnProcessed( [&]( size_t numProcessed, size_t numEnqueued ) {
    progress.Set( ProgressType::WriteHits, numProcessed, numEnqueued );
  } );

  FASTA::Reader< Protein > qryReader( queryPath );

  SequenceList< Protein > queries;
  progress.Activate( ProgressType::ReadQueryFile );
  while( !qryReader.EndOfFile() ) {
    qryReader.Read( numQueriesPerWorkItem, &queries );
    searcher.Enqueue( queries );
    progress.Set( ProgressType::ReadQueryFile, qryReader.NumBytesRead(),
                  qryReader.NumBytesTotal() );
  }

  // Search
  progress.Activate( ProgressType::SearchDB );
  searcher.WaitTillDone();

  progress.Activate( ProgressType::WriteHits );
  writer.WaitTillDone();

  return true;
}
