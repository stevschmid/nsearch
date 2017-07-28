#include "Search.h"

#include <nsearch/Sequence.h>
#include <nsearch/FASTA/Reader.h>
#include <nsearch/Database.h>
#include <nsearch/Database/GlobalSearch.h>

#include "WorkerQueue.h"

#include "Common.h"

using QueryWithHits = std::pair< Sequence, GlobalSearch::HitList >;
using QueryWithHitsList = std::deque< QueryWithHits >;

template<>
class QueueItemInfo< QueryWithHitsList > {
public:
  static size_t Count( const QueryWithHitsList &list ) { return list.size(); }
};

#define MAX_ALIGNMENT_STRING_LENGTH_LINE 50

class SearchResultsWriterWorker {
public:
  SearchResultsWriterWorker( const std::string &path )
    : mWriter( path )
  {

  }

  void Process( const QueryWithHitsList &queryWithHitsList ) {
    for( auto &queryWithHits : queryWithHitsList ) {

      const auto &query = queryWithHits.first;
      const auto &hits = queryWithHits.second;

      // Output with fixed precision (sticky)
      mWriter << std::setiosflags( std::ios::fixed );

      mWriter << "Query >" << query.identifier << std::endl;
      mWriter << " %Id   TLen  Target" << std::endl;
      for( auto &hit : hits ) {
        mWriter
          << std::setprecision( 0 )
          << std::setw( 3 )
          << (hit.alignment.Identity() * 100.0)
          << '%'
          << std::setw( 7 )
          << hit.target.Length()
          << "  >"  << hit.target.identifier
          << std::endl;
      }
      mWriter << std::endl;

      for( auto &hit : hits ) {
        auto queryLen = std::to_string( query.Length() );
        auto targetLen = std::to_string( hit.target.Length() );
        auto maxLen = std::max( queryLen.size(), targetLen.size() );

        mWriter << " Query"
          << std::setw( maxLen + 1 )
          << std::to_string( query.Length() ) << "nt"
          << " >" << query.identifier << std::endl;
        mWriter << "Target"
          << std::setw( maxLen + 1 )
          << std::to_string( hit.target.Length() ) << "nt"
          << " >" << hit.target.identifier << std::endl;


        const Sequence &target = hit.target;
        Cigar cigar = hit.alignment;

        size_t queryStart = 0;
        size_t targetStart = 0;

        // Dont take left terminal gap into account
        if( !cigar.empty() ) {
          const auto &fce = cigar.front();
          if( fce.op == CigarOp::DELETION ) {
            targetStart = fce.count;
            cigar.pop_front();
          } else if( fce.op == CigarOp::INSERTION ) {
            queryStart = fce.count;
            cigar.pop_front();
          }
        }

        // Don't take right terminal gap into account
        if( !cigar.empty() ) {
          const auto &bce = cigar.back();
          if( bce.op == CigarOp::DELETION ) {
            cigar.pop_back();
          } else if( bce.op == CigarOp::INSERTION ) {
            cigar.pop_back();
          }
        }

        bool match;
        size_t numMatches = 0;
        size_t numCols = 0;
        size_t numGaps = 0;

        size_t qcount = queryStart;
        size_t tcount = targetStart;

        bool correct = true;

        using Entry = struct {
          size_t qs, qe;
          std::string q;

          size_t ts, te;
          std::string t;

          std::string a;
        };

        Entry entry;
        entry.qs = 0;
        entry.ts = 0;

        std::deque< Entry > entries;

        for( auto &c : cigar ) {
          for( int i = 0; i < c.count; i++ ) {
            switch( c.op ) {
              case CigarOp::INSERTION:
                entry.t += '-';
                entry.q += query[ qcount++ ];
                entry.a += ' ';
                numGaps++;
                break;

              case CigarOp::DELETION:
                entry.q += '-';
                entry.t += target[ tcount++ ];
                entry.a += ' ';
                numGaps++;
                break;

              case CigarOp::MATCH:
                numMatches++;
                entry.q += query[ qcount++ ];
                entry.t += target[ tcount++ ];
                {
                  bool match = DoNucleotidesMatch( entry.q.back(), entry.t.back() );
                  if( !match ) {
                    correct = false;

                    entry.a += 'X';
                  } else {
                    entry.a += '|';
                  }
                }
                break;

              case CigarOp::MISMATCH:
                entry.a += ' ';
                entry.q += query[ qcount++ ];
                entry.t += target[ tcount++ ];
                break;

              default:
                break;
            }

            numCols++;
            if( numCols % MAX_ALIGNMENT_STRING_LENGTH_LINE == 0 ) {
              entry.qe = qcount;
              entry.te = tcount;
              entries.push_back( entry );

              entry = Entry();
              entry.qs = qcount;
              entry.ts = tcount;
            }
          }
        }

        if( !entry.a.empty() ) {
          entry.qe = qcount;
          entry.te = tcount;
          entries.push_back( entry );
        }

        mWriter << std::endl;
        for( auto &entry : entries ) {
          auto padLen = std::max( std::to_string( queryStart + qcount + 1 ).size(),
              std::to_string( targetStart + tcount + 1 ).size() );

          mWriter << "Qry "
            << std::setw( padLen )
            << entry.qs + queryStart + 1
            << " + " // no strand support for now
            << entry.q
            << " "
            << entry.qe + queryStart
            << std::endl;

          mWriter << std::string( 9, ' ' )
            << entry.a
            << std::endl;

          mWriter << "Tgt "
            << std::setw( padLen )
            << entry.ts + targetStart + 1
            << " + " // no strand support for now
            << entry.t
            << " "
            << entry.te + targetStart
            << std::endl;

          mWriter << std::endl;
        }

        if( !correct ) {
          mWriter << "[DEBUG] INVALID ALIGNMENT" << std::endl;
        }

        float identity = float( numMatches ) / float( numCols );
        float gapsRatio = float( numGaps ) / float( numCols );
        mWriter << numCols << " cols, "
          << numMatches << " ids ("
          << std::setprecision( 1 ) << ( 100.0f * identity ) << "%), "
          << numGaps << " gaps ("
          << std::setprecision( 1 ) << ( 100.0f * gapsRatio ) << "%)" << std::endl;

        mWriter << std::endl;
      }

    }
  }

private:
  std::ofstream mWriter;
};
using SearchResultsWriter = WorkerQueue< SearchResultsWriterWorker, QueryWithHitsList, const std::string& >;

template<>
class QueueItemInfo< SequenceList > {
public:
  static size_t Count( const SequenceList &list ) { return list.size(); }
};

class QueryDatabaseSearcherWorker {
public:
  QueryDatabaseSearcherWorker( SearchResultsWriter &writer, const Database &database )
    : mWriter( writer), mGlobalSearch( database, 0.75, 1, 8 )
  {
  }

  void Process( const SequenceList &queries ) {
    QueryWithHitsList list;

    for( auto &query : queries ) {
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
  GlobalSearch mGlobalSearch;
  SearchResultsWriter &mWriter;
};
using QueryDatabaseSearcher = WorkerQueue< QueryDatabaseSearcherWorker, SequenceList, SearchResultsWriter&, const Database& >;

bool Search( const std::string &queryPath, const std::string &databasePath, const std::string &outputPath ) {
  ProgressOutput progress;

  Sequence seq;
  SequenceList sequences;

  FASTA::Reader dbReader( databasePath );

  enum ProgressType { ReadDBFile, StatsDB, IndexDB, ReadQueryFile, SearchDB, WriteHits };

  progress.Add( ProgressType::ReadDBFile, "Reading DB file", UnitType::BYTES );
  progress.Add( ProgressType::StatsDB, "Analyzing DB sequences");
  progress.Add( ProgressType::IndexDB, "Indexing database");
  progress.Add( ProgressType::ReadQueryFile, "Reading query file", UnitType::BYTES );
  progress.Add( ProgressType::SearchDB, "Searching database" );
  progress.Add( ProgressType::WriteHits, "Write hits" );

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

  // Read and process queries
  const int numQueriesPerWorkItem = 50;

  SearchResultsWriter writer( 1, outputPath );
  QueryDatabaseSearcher searcher( -1, writer, db );

  searcher.OnProcessed( [&]( size_t numProcessed, size_t numEnqueued ) {
    progress.Set( ProgressType::SearchDB, numProcessed, numEnqueued );
  });

  writer.OnProcessed( [&]( size_t numProcessed, size_t numEnqueued ) {
    progress.Set( ProgressType::WriteHits, numProcessed, numEnqueued );
  });

  FASTA::Reader qryReader( queryPath );

  SequenceList queries;
  progress.Activate( ProgressType::ReadQueryFile );
  while( !qryReader.EndOfFile() )  {
    qryReader.Read( queries, numQueriesPerWorkItem );
    searcher.Enqueue( queries );
    progress.Set( ProgressType::ReadQueryFile, qryReader.NumBytesRead(), qryReader.NumBytesTotal() );
  }

  // Search
  progress.Activate( ProgressType::SearchDB );
  searcher.WaitTillDone();

  progress.Activate( ProgressType::WriteHits );
  writer.WaitTillDone();

  return true;
}
