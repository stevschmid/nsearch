#pragma once

#include <sparsepp/spp.h>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <cassert>
#include <iostream>
#include <iomanip>
#include <memory>

#include "Sequence.h"
#include "Utils.h"
#include "Aligner.h"
#include "SpacedSeedsSIMD.h"

#include "Alignment/Seed.h"
#include "Alignment/Ranges.h"
#include "Alignment/HitTracker.h"
#include "Alignment/ExtendAlign.h"
#include "Alignment/BandedAlign.h"

// HSP: first and last character in sequence (i.e. seq[a1] - seq[a2])
class HSP {
public:
  size_t a1, a2;
  size_t b1, b2;
  int score;
  Cigar cigar;

  HSP( size_t a1, size_t a2, size_t b1, size_t b2 )
    : a1( a1 ), a2( a2 ), b1( b1 ), b2( b2 ), score( 0 )
  {
    assert( a2 >= a1 && b2 >= b1 );
  }

  size_t Length() const {
    return std::max( a2 - a1, b2 - b1 ) + 1;
  }

  int Score() const {
    return score;
  }

  bool IsOverlapping( const HSP &other ) const {
    return
      ( a1 <= other.a2 && other.a1 <= a2 ) // overlap in A direction
      ||
      ( b1 <= other.b2 && other.b1 <= b2 ); // overlap in B direction
  }

  size_t DistanceTo( const HSP &other ) const {
    size_t x = ( a1 > other.a2 ? a1 - other.a2 : other.a1 - a2 );
    size_t y = ( b1 > other.b2 ? b1 - other.b2 : other.b2 - b2 );

    return sqrt( x*x + y*y );
  }

  bool operator<( const HSP &other ) const {
    return Score() < other.Score();
  }
};

bool PrintWholeAlignment( const Sequence &query, const Sequence &target, const Cigar &cigar_ ) {
  std::string q;
  std::string t;
  std::string a;
  bool correct = true;

  size_t queryStart = 0;
  size_t targetStart = 0;

  Cigar cigar = cigar_;

  // Dont take left terminal gap into account
  if( !cigar.empty() ) {
    const CigarEntry &fce = cigar.front();
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
    const CigarEntry &bce = cigar.back();
    if( bce.op == CigarOp::DELETION ) {
      cigar.pop_back();
    } else if( bce.op == CigarOp::INSERTION ) {
      cigar.pop_back();
    }
  }

  bool match;
  size_t numMatches = 0;
  size_t numCols = 0;

  size_t qcount = queryStart;
  size_t tcount = targetStart;

  for( auto &c : cigar ) {
    for( int i = 0; i < c.count; i++ ) {
      switch( c.op ) {
        case CigarOp::INSERTION:
          t += '-';
          q += query[ qcount++ ];
          a += ' ';
          break;

        case CigarOp::DELETION:
          q += '-';
          t += target[ tcount++ ];
          a += ' ';
          break;

        case CigarOp::MATCH:
          numMatches++;
          q += query[ qcount++ ];
          t += target[ tcount++ ];
          {
            bool match = DoNucleotidesMatch( q.back(), t.back() );
            if( !match ) {
              correct = false;
              a += 'X';
            } else {
              a += '|';
            }
          }
          break;

        case CigarOp::MISMATCH:
          a += ' ';
          q += query[ qcount++ ];
          t += target[ tcount++ ];
          break;

        default:
          break;
      }

      numCols++;
    }
  }

  std::cout << std::endl;
  std::cout << "Query " << std::string( 11, ' ' ) << ">" << query.identifier << std::endl;
  std::cout << std::endl;
  std::cout << std::setw( 15 ) << queryStart + 1 << " " << q << " " << qcount << std::endl;
  std::cout << std::string( 16, ' ' ) << a << std::endl;
  std::cout << std::setw( 15 ) << targetStart + 1 << " " << t << " " << tcount << std::endl;
  std::cout << std::endl;
  std::cout << "Target " << std::string( 11, ' ' ) << ">" << target.identifier << std::endl;

  std::cout << std::endl;
  float identity = float( numMatches ) / float( numCols );
  std::cout <<  numCols << " cols, " << numMatches << " ids (" << ( 100.0f * identity ) << "%)" << std::endl;

  std::cout << std::endl;
  std::cout << std::string( 50, '=') << std::endl;

  return correct;
}

class Database {

  float CalculateIdentity( const Cigar &cigar ) const {
    size_t cols = 0;
    size_t matches = 0;

    for( const CigarEntry &c : cigar ) {
      // Don't count terminal gaps towards identity calculation
      if( &c == &(*cigar.cbegin()) && ( c.op == CigarOp::INSERTION || c.op == CigarOp::DELETION ) )
        continue;
      if( &c == &(*cigar.crbegin()) && ( c.op == CigarOp::INSERTION || c.op == CigarOp::DELETION ) )
        continue;

      cols += c.count;
      if( c.op == CigarOp::MATCH )
        matches += c.count;
    }

    return cols > 0 ? float( matches ) / float( cols ) : 0.0f;
  }

public:

  void Stats() const {
  }

  Database( const SequenceList &sequences, size_t wordSize )
    : mSequences( sequences ), mWordSize( wordSize )
  {
    mNumUniqueWords = 1 << ( 2 * mWordSize ); // 2 bits per nt

    mTotalNucleotides = 0;

    mTotalWords = 0;
    for( auto &seq : mSequences ) {
      mTotalWords += seq.Length() - mWordSize + 1;
    }

    using Entry = struct Entry_s {
      uint32_t word;
      uint32_t index;
      uint32_t pos;

      Entry_s(uint32_t word, uint32_t index, uint32_t pos )
        : word( word ), index( index ), pos( pos )
      {
      }
    };

    std::vector< Entry > entries;
    entries.reserve( mTotalWords );

    size_t totalEntries = 0;
    size_t totalFirstEntries = 0;
    std::vector< uint32_t > uniqueCount( mNumUniqueWords );
    std::vector< uint32_t > uniqueIndex( mNumUniqueWords, -1 );
    for( uint32_t idx = 0; idx < mSequences.size(); idx++ ) {
      const Sequence &seq = mSequences[ idx ];
      mTotalNucleotides += seq.Length();

      SpacedSeedsSIMD spacedSeeds( seq, mWordSize );
      spacedSeeds.ForEach( [&]( size_t pos, uint32_t word ) {
        entries.emplace_back( word, idx, pos );
        totalEntries++;

        if( uniqueIndex[ word ] != idx ) {
          uniqueIndex[ word ] = idx;
          uniqueCount[ word ]++;
          totalFirstEntries++;
        }
      });
    }
    mTotalWords = entries.size();

    // Calculate indices
    /* mWordIndices.reserve( mNumUniqueWords ); */
    mIndexByWord.reserve( mNumUniqueWords );
    for( size_t i = 0; i < mNumUniqueWords; i++ ) {
      mIndexByWord[ i ] = i > 0 ? mIndexByWord[ i - 1 ] + uniqueCount[ i - 1 ] : 0;
    }

    // Populate DB
    mFirstEntries.resize( totalFirstEntries );
    mFurtherEntries.reserve( totalEntries - totalFirstEntries );

    // Entries is sorted by sequence
    mNumEntriesByWord = std::vector< uint32_t >( mNumUniqueWords );
    for( auto &s : entries ) {
      uint32_t word = s.word;

      if( mNumEntriesByWord[ word ] == 0 ) {
        // Create new entry
        WordEntry *entry = &mFirstEntries[ mIndexByWord[ word ] ];
        entry->sequence = s.index;
        entry->pos = s.pos;
        entry->nextEntry = NULL;
        mNumEntriesByWord[ word ]++;
      } else {
        // Check if last entry == index
        WordEntry *entry = &mFirstEntries[ mIndexByWord[ word ] + ( mNumEntriesByWord[ word ] - 1 ) ];

        if( entry->sequence == s.index ) {
          mFurtherEntries.emplace_back( s.index, s.pos );
          WordEntry *s = entry;
          while( s->nextEntry ) {
            s = s->nextEntry;
          }
          s->nextEntry = &mFurtherEntries.back();
        } else {
          entry++;
          entry->sequence = s.index;
          entry->pos = s.pos;
          entry->nextEntry = NULL;
          mNumEntriesByWord[ word ]++;
        }
      }
    }
  }

  SequenceList Query( const Sequence &query, float minIdentity, int maxHits = 1, int maxRejects = 8 ) {
    const size_t defaultMinHSPLength = 16;
    const size_t maxHSPJoinDistance = 16;

    std::cout << "===> SEARCH " << query.identifier << std::endl;

    size_t minHSPLength = std::min( defaultMinHSPLength, query.Length() / 2 );

    // Go through each kmer, find hits
    if( mHits.size() < mSequences.size() ) {
      mHits.resize( mSequences.size() );
    }

    // Fast counter reset
    memset( mHits.data(), 0, sizeof( size_t ) * mHits.capacity() );

    class Candidate {
    public:
      size_t seqIdx;
      size_t counter;

      Candidate( size_t seqIdx, size_t counter )
        : seqIdx( seqIdx ), counter( counter )
      {
      }

      bool operator<( const Candidate &other ) const {
        return counter < other.counter;
      }
    };

    // TODO: Fix mWorDSize in AddHIT!
    // TODO: Cache SPACED SEEDS! similar to db()

    std::multiset< Candidate > highscores;
    SpacedSeedsSIMD spacedSeeds( query, mWordSize );
    std::vector< bool > uniqueCheck( mNumUniqueWords );

    using Kmer = struct Kmer_s {
      size_t pos;
      uint32_t word;

      Kmer_s( size_t p, uint32_t w )
        : pos( p ), word( w )
      {
      }
    };
    std::vector< Kmer > kmers;
    kmers.reserve( query.Length() );

    spacedSeeds.ForEach( [&]( size_t pos, uint32_t word ) {
      kmers.emplace_back( pos, word );

      if( !uniqueCheck[ word ] ) {
        uniqueCheck[ word ] = 1;

        WordEntry *ptr = &mFirstEntries[ mIndexByWord[ word ] ];
        for( uint32_t i = 0; i < mNumEntriesByWord[ word ]; i++, ptr++ ) {
          uint32_t candidateIdx = ptr->sequence;

          size_t &counter = mHits[ candidateIdx ];
          counter++;

          if( highscores.empty() || counter > highscores.begin()->counter ) {
            auto it = std::find_if( highscores.begin(), highscores.end(), [candidateIdx]( const Candidate &c ) {
              return c.seqIdx == candidateIdx;
            });

            if( it != highscores.end() ) {
              highscores.erase( it );
            }

            highscores.insert( Candidate( candidateIdx, counter ) );
            if( highscores.size() > maxHits + maxRejects ) {
              highscores.erase( highscores.begin() );
            }
          }
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
    std::cout << "Highscores " << highscores.size() << std::endl;
    for( auto it = highscores.rbegin(); it != highscores.rend(); ++it ) {
      const size_t seqIdx = (*it).seqIdx;
      assert( seqIdx < mSequences.count() );
      const Sequence &candidateSeq = mSequences[ seqIdx ];
      std::cout << "Highscore Entry " << seqIdx << " " << (*it).counter << std::endl;

      // Go through each kmer, find hits
      HitTracker hitTracker;

      for( auto &k : kmers ) {
        WordEntry *ptr = &mFirstEntries[ mIndexByWord[ k.word ] ];

        for( uint32_t i = 0; i < mNumEntriesByWord[ k.word ]; i++, ptr++ ) {
          if( ptr->sequence != seqIdx )
            continue;

          while( ptr ) {
            hitTracker.AddHit( k.pos, ptr->pos, spacedSeeds.Length() );
            ptr = ptr->nextEntry;
          }
          break;
        }
      }

      // Find all HSP
      // Sort by length
      // Try to find best chain
      // Fill space between with banded align

      std::set< HSP > hsps;
      for( auto &seed : hitTracker.Seeds() ) {
        size_t queryPos, candidatePos;


        size_t a1 = seed.s1, a2 = seed.s1 + seed.length - 1,
               b1 = seed.s2, b2 = seed.s2 + seed.length - 1;

        Cigar leftCigar;
        int leftScore = mExtendAlign.Extend( query, candidateSeq,
            &queryPos, &candidatePos,
            &leftCigar,
            AlignmentDirection::backwards,
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
            AlignmentDirection::forwards,
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
          for( size_t s1 = seed.s1, s2 = seed.s2;
               s1 < seed.s1 + seed.length && s2 < seed.s2 + seed.length;
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
        mBandedAlign.Align( query, candidateSeq, &cigar, AlignmentDirection::backwards, first.a1, first.b1 );
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
              AlignmentDirection::forwards,
              current.a2 + 1, current.b2 + 1,
              next.a1, next.b1 );
          alignment += cigar;
        }

        // Align last HSP's end to whole sequences end
        auto &last = *chain.crbegin();
        alignment += last.cigar;
        mBandedAlign.Align( query, candidateSeq, &cigar, AlignmentDirection::forwards, last.a2 + 1, last.b2 + 1 );
        alignment += cigar;

        float identity = CalculateIdentity( alignment );
        if( identity >= minIdentity ) {
          accept = true;
          bool correct = PrintWholeAlignment( query, candidateSeq, alignment );
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

    return SequenceList();
  }

private:
  ExtendAlign mExtendAlign;
  BandedAlign mBandedAlign;

  size_t mWordSize;

  std::vector< size_t > mHits;

  SequenceList mSequences;
  size_t mTotalNucleotides;
  size_t mTotalWords;
  size_t mNumUniqueWords;

  /* using WordInfo = struct { */
  /*   uint32_t index; */
  /*   uint32_t pos; */
  /* }; */
  /* using WordInfoDatabase = std::vector< WordInfo >; */

  /* std::vector< uint32_t > mWordIndices; */
  /* std::vector< uint32_t > mWordCounts; */
  /* WordInfoDatabase mWordInfoDB; */

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
