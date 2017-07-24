#pragma once

#include <unordered_map>
#include <unordered_set>
#include <set>
#include <cassert>
#include <iostream>
#include <iomanip>
#include <memory>

#include "Sequence.h"
#include "Utils.h"

#include "Alignment/Common.h"
#include "Alignment/Ranges.h"
#include "Alignment/HitTracker.h"
#include "Alignment/ExtendAlign.h"
#include "Alignment/BandedAlign.h"

#include "Database/Kmers.h"
#include "Database/HSP.h"

class Highscore {
  class Entry {
  public:
    size_t id = 0;
    size_t score = 0;

    bool operator<( const Entry &other ) const {
      return score < other.score;
    }
  };
public:
  Highscore( size_t numHighestEntriesToKeep )
    : mLowestScore( 0 )
  {
    mEntries.resize( numHighestEntriesToKeep );
  }

  // score is assumed to increase for every id
  void Set( size_t id, size_t score ) {
    if( score < mLowestScore )
      return;

    auto it = std::find_if( mEntries.begin(), mEntries.end(), [ id ]( const Entry &candidate ) {
        return id == candidate.id;
        });

    if( it == mEntries.end() ) {
      it = std::find_if( mEntries.begin(), mEntries.end(), [ score ]( const Entry &candidate ) {
          return score > candidate.score;
          });
    }

    if( it != mEntries.end() ) {
      it->id = id;
      it->score = score;

      mLowestScore = std::min_element( mEntries.begin(), mEntries.end() )->score;
    }
  }

  std::vector< Entry > EntriesFromTopToBottom() const {
    std::vector< Entry > topToBottom = mEntries;
    std::sort( topToBottom.begin(), topToBottom.end(), []( const Entry &a, const Entry &b ) {
        return !(a < b);
        });
    return topToBottom;
  }

private:
  size_t mLowestScore;
  std::vector< Entry > mEntries;
};

class DatabaseSearcher;

class Database {
  friend class DatabaseSearcher;

public:

  size_t Size() const {
    return mSequences.size();
  }

  void Stats() const {
  }

  Database( const SequenceList &sequences, size_t wordSize )
    : mSequences( sequences ), mWordSize( wordSize )
  {
    mMaxUniqueWords = 1 << ( 2 * mWordSize ); // 2 bits per nt

    size_t totalEntries = 0;
    size_t totalFirstEntries = 0;
    std::vector< uint32_t > uniqueCount( mMaxUniqueWords );
    std::vector< uint32_t > uniqueIndex( mMaxUniqueWords, -1 );
    for( uint32_t idx = 0; idx < mSequences.size(); idx++ ) {
      const Sequence &seq = mSequences[ idx ];

      Kmers spacedSeeds( seq, mWordSize );
      spacedSeeds.ForEach( [&]( Kmer word, size_t pos ) {
        totalEntries++;

        if( uniqueIndex[ word ] != idx ) {
          uniqueIndex[ word ] = idx;
          uniqueCount[ word ]++;
          totalFirstEntries++;
        }
      });
    }

    // Calculate indices
    mIndexByWord.reserve( mMaxUniqueWords );
    for( size_t i = 0; i < mMaxUniqueWords; i++ ) {
      mIndexByWord[ i ] = i > 0 ? mIndexByWord[ i - 1 ] + uniqueCount[ i - 1 ] : 0;
    }

    // Populate DB
    mFirstEntries.resize( totalFirstEntries );
    mFurtherEntries.reserve( totalEntries - totalFirstEntries );

    // Entries is sorted by sequence
    mNumEntriesByWord = std::vector< uint32_t >( mMaxUniqueWords );
    for( uint32_t idx = 0; idx < mSequences.size(); idx++ ) {
      const Sequence &seq = mSequences[ idx ];

      Kmers spacedSeeds( seq, mWordSize );
      spacedSeeds.ForEach( [&]( Kmer word, size_t pos ) {
        if( mNumEntriesByWord[ word ] == 0 ) {
          // Create new entry
          WordEntry *entry = &mFirstEntries[ mIndexByWord[ word ] ];
          entry->sequence = idx;
          entry->pos = pos;
          entry->nextEntry = NULL;
          mNumEntriesByWord[ word ]++;
        } else {
          // Check if last entry == index
          WordEntry *entry = &mFirstEntries[ mIndexByWord[ word ] + ( mNumEntriesByWord[ word ] - 1 ) ];

          if( entry->sequence == idx ) {
            mFurtherEntries.emplace_back( idx, pos );
            WordEntry *s = entry;
            while( s->nextEntry ) {
              s = s->nextEntry;
            }
            s->nextEntry = &mFurtherEntries.back();
          } else {
            entry++;
            entry->sequence = idx;
            entry->pos = pos;
            entry->nextEntry = NULL;
            mNumEntriesByWord[ word ]++;
          }
        }
      });

    }
  }

private:
  ExtendAlign mExtendAlign;
  BandedAlign mBandedAlign;

  size_t mWordSize;

  std::vector< size_t > mHits;

  SequenceList mSequences;
  size_t mMaxUniqueWords;

  using WordEntry = struct WordEntry_s {
    uint32_t sequence;
    uint32_t pos;
    WordEntry_s *nextEntry;

    WordEntry_s()
    : sequence( -1 )
    {
    }

    WordEntry_s( uint32_t s, uint32_t p, WordEntry_s *ne = NULL )
      : sequence( s ), pos( p ), nextEntry( ne )
    {

    }
  };

  std::vector< uint32_t > mIndexByWord;
  std::vector< uint32_t > mNumEntriesByWord;
  std::vector< WordEntry > mFirstEntries; // first word (kmer) hit for each candidate
  std::vector< WordEntry > mFurtherEntries;
};

class DatabaseSearcher {
public:
  using Result = struct {
    Sequence query;
    Sequence target;

    Cigar alignment;

    size_t targetStart, targetLength;
    size_t queryStart, queryLength;

    size_t numCols, numMatches, numMismatches, numGaps;

    float identity;
  };
  using ResultList = std::deque< Result >;

  DatabaseSearcher( const Database &db )
    : mDB( db )
  {
  }

  ResultList Query( const Sequence &query, float minIdentity, int maxHits = 1, int maxRejects = 8 ) {

    const size_t defaultMinHSPLength = 16;
    const size_t maxHSPJoinDistance = 16;

    std::cout << "===> SEARCH " << query.identifier << std::endl;

    size_t minHSPLength = std::min( defaultMinHSPLength, query.Length() / 2 );

    // Go through each kmer, find hits
    if( mHits.size() < mDB.Size() ) {
      mHits.resize( mDB.Size() );
    }

    // Fast counter reset
    memset( mHits.data(), 0, sizeof( size_t ) * mHits.capacity() );

    Highscore highscore( maxHits + maxRejects );

    Kmers spacedSeeds( query, mDB.mWordSize );
    std::vector< bool > uniqueCheck( mDB.mMaxUniqueWords );

    spacedSeeds.ForEach( [&]( Kmer word, size_t pos ) {
      if( !uniqueCheck[ word ] ) {
        uniqueCheck[ word ] = 1;

        const Database::WordEntry *ptr = &mDB.mFirstEntries[ mDB.mIndexByWord[ word ] ];
        for( uint32_t i = 0; i < mDB.mNumEntriesByWord[ word ]; i++, ptr++ ) {
          uint32_t candidateIdx = ptr->sequence;

          size_t &counter = mHits[ candidateIdx ];
          counter++;

          highscore.Set( candidateIdx, counter );
        }
      }
    });

    // For each candidate:
    // - Get HSPs,
    // - Check for good HSP (>= similarity threshold)
    // - Join HSP together
    // - Align
    // - Check similarity
    int numHits = 0;
    int numRejects = 0;

    auto highscores = highscore.EntriesFromTopToBottom();
    std::cout << "Highscores " << highscores.size() << std::endl;

    for( auto it = highscores.cbegin(); it != highscores.cend(); ++it ) {
      const size_t seqIdx = it->id;
      assert( seqIdx < mDB.mSequences.size() );
      const Sequence &candidateSeq = mDB.mSequences[ seqIdx ];
      std::cout << "Highscore Entry " << it->id << " " << it->score << std::endl;

      // Go through each kmer, find hits
      HitTracker hitTracker;

      spacedSeeds.ForEach( [&]( Kmer word, size_t pos ) {
        const Database::WordEntry *ptr = &mDB.mFirstEntries[ mDB.mIndexByWord[ word ] ];
        for( uint32_t i = 0; i < mDB.mNumEntriesByWord[ word ]; i++, ptr++ ) {
          if( ptr->sequence != seqIdx )
            continue;

          while( ptr ) {
            hitTracker.AddHit( pos, ptr->pos, spacedSeeds.Length() );
            ptr = ptr->nextEntry;
          }
          break;
        }
      });

      // Find all HSP
      // Sort by length
      // Try to find best chain
      // Fill space between with banded align

      std::set< HSP > hsps;
      for( auto &sp : hitTracker.List() ) {
        size_t queryPos, candidatePos;

        size_t a1 = sp.s1, a2 = sp.s1 + sp.length - 1,
               b1 = sp.s2, b2 = sp.s2 + sp.length - 1;

        Cigar leftCigar;
        int leftScore = mExtendAlign.Extend( query, candidateSeq,
            &queryPos, &candidatePos,
            &leftCigar,
            AlignmentDirection::rev,
            a1, b1 );
        if( !leftCigar.empty() ) {
          a1 = queryPos;
          b1 = candidatePos;
        }

        Cigar rightCigar;
        size_t rightQuery, rightCandidate;
        int rightScore = mExtendAlign.Extend( query, candidateSeq,
            &queryPos, &candidatePos,
            &rightCigar,
            AlignmentDirection::fwd,
            a2 + 1, b2 + 1 );
        if( !rightCigar.empty() ) {
          a2 = queryPos;
          b2 = candidatePos;
        }

        HSP hsp( a1, a2, b1, b2 );
        if( hsp.Length() >= minHSPLength ) {
          // Construct hsp cigar (spaced seeds so we cannot assume full match)
          Cigar middleCigar;
          int middleScore = 0;
          for( size_t s1 = sp.s1, s2 = sp.s2;
               s1 < sp.s1 + sp.length && s2 < sp.s2 + sp.length;
               s1++, s2++ )
          {
            bool match = DoNucleotidesMatch( query[ s1 ], candidateSeq[ s2 ] );
            middleCigar.Add( match ? CigarOp::MATCH : CigarOp::MISMATCH );
            middleScore += match ? mExtendAlign.AP().matchScore : mExtendAlign.AP().mismatchScore;
          }
          hsp.score = leftScore + middleScore + rightScore;
          hsp.cigar = leftCigar + middleCigar + rightCigar;

          // Save HSP
          hsps.insert( hsp );
        }
      }

      // Greedy join HSPs if close
      struct HSPChainOrdering {
        bool operator() ( const HSP &left, const HSP &right ) const {
          return left.a1 < right.a1 && left.b1 < right.b1;
        }
      };

      std::set< HSP, HSPChainOrdering > chain;
      for( auto it = hsps.rbegin(); it != hsps.rend(); ++it ) {
        const HSP &hsp = *it;
        bool hasNoOverlaps = std::none_of( chain.begin(), chain.end(), [&]( const HSP &existing ) {
          return hsp.IsOverlapping( existing );
        });
        if( hasNoOverlaps ) {
          bool anyHSPJoinable = std::any_of( chain.begin(), chain.end(), [&]( const HSP &existing ) {
            return hsp.DistanceTo( existing ) <= maxHSPJoinDistance;
          });

          if( chain.empty() || anyHSPJoinable ) {
            chain.insert( hsp );
          }
        }
      }

      bool accept = false;
      if( chain.size() > 0 ) {
        Cigar alignment;
        Cigar cigar;

        // Align first HSP's start to whole sequences begin
        auto &first = *chain.cbegin();
        mBandedAlign.Align( query, candidateSeq, &cigar, AlignmentDirection::rev, first.a1, first.b1 );
        alignment += cigar;

        // Align in between the HSP's
        for( auto it1 = chain.cbegin(), it2 = ++chain.cbegin();
            it1 != chain.cend() && it2 != chain.cend();
            ++it1, ++it2 )
        {
          auto &current = *it1;
          auto &next = *it2;

          alignment += current.cigar;
          mBandedAlign.Align( query, candidateSeq, &cigar,
              AlignmentDirection::fwd,
              current.a2 + 1, current.b2 + 1,
              next.a1, next.b1 );
          alignment += cigar;
        }

        // Align last HSP's end to whole sequences end
        auto &last = *chain.crbegin();
        alignment += last.cigar;
        mBandedAlign.Align( query, candidateSeq, &cigar, AlignmentDirection::fwd, last.a2 + 1, last.b2 + 1 );
        alignment += cigar;

        float identity = alignment.Identity();
        if( identity >= minIdentity ) {
          accept = true;

          bool correct;

          std::cout << std::endl;
          std::cout << alignment.ToFullAlignmentString( query, candidateSeq, &correct );
          std::cout << std::endl;
          std::cout << std::string( 50, '=') << std::endl;

          if( !correct ) {
            std::cout << "INVALID ALIGNMENT" << std::endl;
          }
        }
      }

      if( accept ) {
        numHits++;
        if( numHits >= maxHits )
          break;
      } else {
        numRejects++;
        if( numRejects >= maxRejects )
          break;
      }
    }

    return ResultList();
  }

private:
  const Database &mDB;

  std::vector< size_t > mHits;

  ExtendAlign mExtendAlign;
  BandedAlign mBandedAlign;
};
