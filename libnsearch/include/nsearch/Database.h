#pragma once

#include <unordered_map>
#include <unordered_set>
#include <set>

#include "Sequence.h"
#include "Utils.h"
#include "Aligner.h"
#include "Kmer.h"

#include "Alignment/Seed.h"
#include "Alignment/Ranges.h"
#include "Alignment/HitTracker.h"
#include "Alignment/ExtendAlign.h"
#include "Alignment/BandedAlign.h"
#include "Alignment/OptimalChainFinder.h"

class Database {
  using SequenceRef = std::shared_ptr< Sequence >;
  using SequenceInfo = std::pair< size_t, SequenceRef >;
  using SequenceMappingDatabase = std::unordered_map< Sequence, std::deque< SequenceInfo > >;

public:
  Database( size_t wordSize )
    : mWordSize( wordSize )
  {
    // Sequences have to fit in the calculated hash
    assert( mWordSize * 2 <= sizeof( size_t ) * 8 );
  }

  void AddSequence( const Sequence &seq ) {
    // Save
    mSequences.push_back( std::make_shared< Sequence >( seq ) );
    SequenceRef ref = mSequences.back();

    // Kmers for Indexing
    Kmers kmers( *ref, mWordSize );
    kmers.ForEach( [&]( const Sequence &kmer, size_t pos ) {
      this->mWordDB[ kmer ].push_back( SequenceInfo( pos, ref ) );
    });
  }

  SequenceList Query( const Sequence &query, int maxHits = 10 ) {
    const int xDrop = 32;
    const size_t defaultMinHSPLength = 16;
    const size_t maxHSPJoinDistance = 16;

    size_t minHSPLength = std::min( defaultMinHSPLength, query.Length() / 2 );

    std::unordered_map< SequenceRef, HitTracker > hits;

    // Go through each kmer, find hits
    Kmers kmers( query, mWordSize );
    kmers.ForEach( [&]( const Sequence &kmer, size_t pos ) {
      for( auto &seqInfo : mWordDB[ kmer ] ) {
        size_t candidatePos = seqInfo.first;
        SequenceRef candidateRef = seqInfo.second;
        hits[ candidateRef ].AddHit( pos, candidatePos, mWordSize );
      }
    });

    // Order candidates by number of shared words (hits)
    using Candidate = std::pair< SequenceRef, HitTracker >;
    struct CandidateCompare {
      bool operator() ( const Candidate &left, const Candidate &right ) const {
        return left.second.Score() < right.second.Score();
      }
    };
    using Highscore = std::set< Candidate , CandidateCompare >;
    Highscore highscore;
    for( auto &h : hits ) {
      highscore.insert( std::make_pair( h.first, h.second ) );
    }

    // For each candidate:
    // - Get HSPs,
    // - Check for good HSP (>= similarity threshold)
    // - Join HSP together
    // - Align
    // - Check similarity
    /* std::cout << "=====NEWQUERY=====" << std::endl; */
    /* std::cout << "MinHSP Length " << minHSPLength << std::endl; */

    for( auto it = highscore.rbegin(); it != highscore.rend(); ++it ) {
      const Candidate &candidate = *it;
      const Sequence &candidateSeq = *candidate.first;

      // some yolo preliminary optimization
      /* if( candidate.second.Score() < 0.5 * query.Length() ) { */
      /*   continue; */
      /* } */

      // Find all HSP
      // Sort by length
      // Try to find best chain
      // Fill space between with banded align

      std::multimap< size_t, Seed > hsps;

      for( auto &seed : candidate.second.Seeds() ) {
        size_t leftQuery, leftCandidate;
        int leftScore = mExtendAlign.Extend( query, candidateSeq,
            xDrop,
            AlignmentDirection::backwards,
            seed.s1, seed.s2,
            &leftQuery, &leftCandidate );

        size_t rightQuery, rightCandidate;
        int rightScore = mExtendAlign.Extend( query, candidateSeq,
            xDrop,
            AlignmentDirection::forwards,
            seed.s1 + seed.length, seed.s2 + seed.length,
            &rightQuery, &rightCandidate );

        Seed hsp = seed;
        hsp.s1 = leftQuery;
        hsp.s2 = leftCandidate;
        hsp.length = rightQuery - leftQuery;

        if( hsp.length >= minHSPLength ) {
          hsps.insert( std::make_pair( hsp.length, hsp ) );
        }
      }

      // Greedy join HSPs if close
      auto precede = []( const Seed &seed, const Seed &other ) {
        return seed.s1 + seed.length <= other.s1 &&
          seed.s2 + seed.length <= other.s2;
      };

      auto succeed = []( const Seed &seed, const Seed &other ) {
        return seed.s1 >= other.s1 + other.length &&
          seed.s2 >= other.s2 + other.length;
      };

      struct SeedCompare {
        bool operator() ( const Seed &left, const Seed &right ) const {
          return left.s1 < right.s1 && left.s2 < right.s2;
        }
      };

      std::set< Seed, SeedCompare > chain;
      for( auto &p : hsps ) {
        const Seed &hsp = p.second;
        if( chain.empty() ) {
          chain.insert( hsp );
        } else {
          const Seed &left = *chain.begin();
          const Seed &right = *chain.rbegin();

          if( precede( hsp, left ) ) {
            size_t dist = std::max( left.s1 - ( hsp.s1 + hsp.length ),
                left.s2 - ( hsp.s2 + hsp.length ) );

            if( dist <= maxHSPJoinDistance ) {
              chain.insert( hsp );
            }
          } else if( succeed( hsp, right ) ) {
            size_t dist = std::max( hsp.s1 - ( right.s1 + right.length ),
                hsp.s2 - ( right.s2 + right.length ) );

            if( dist <= maxHSPJoinDistance ) {
              chain.insert( hsp );
            }
          }
        }
      }

      if( chain.size() > 0 ) {
        Cigar alignment;
        Cigar cigar;

        // Align first HSP's start to whole sequences begin
        auto &first = *chain.cbegin();
        mBandedAlign.Align( query, candidateSeq, &cigar, first.s1, first.s2, AlignmentDirection::backwards );
        alignment += cigar;

        // Align in between the HSP's
        for( auto it1 = chain.cbegin(), it2 = ++chain.cbegin();
            it1 != chain.cend() && it2 != chain.cend();
            ++it1, ++it2 )
        {
          auto &prev = *it1;
          auto &next = *it2;

          mBandedAlign.Align( query, candidateSeq, &cigar,
              prev.s1 + prev.length, prev.s2 + prev.length,
              AlignmentDirection::forwards,
              next.s1, next.s2 );
          alignment += cigar;
        }

        // Align last HSP's end to whole sequences end
        auto &last = *chain.crbegin();
        mBandedAlign.Align( query, candidateSeq, &cigar, last.s1 + last.length, last.s2 + last.length, AlignmentDirection::forwards );
        alignment += cigar;
      }
    }
    return SequenceList();
  }

private:
  ExtendAlign mExtendAlign;
  BandedAlign mBandedAlign;

  size_t mWordSize;
  std::deque< SequenceRef > mSequences;

  SequenceMappingDatabase mWordDB;
};
