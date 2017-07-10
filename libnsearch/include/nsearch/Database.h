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
#include "Aligner.h"
#include "HashWords.h"

#include "Alignment/Seed.h"
#include "Alignment/Ranges.h"
#include "Alignment/HitTracker.h"
#include "Alignment/ExtendAlign.h"
#include "Alignment/BandedAlign.h"
#include "Alignment/OptimalChainFinder.h"

class HSP {
public:
  size_t a1, a2;
  size_t b1, b2;
  Cigar cigar;

  HSP( size_t a1, size_t a2, size_t b1, size_t b2 )
    : a1( a1 ), a2( a2 ), b1( b1 ), b2( b2 )
  {
    assert( a2 >= a1 && b2 >= b1 );
  }

  size_t Length() const {
    return std::max( a2 - a1, b2 - b1 );
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
    return Length() < other.Length();
  }
};

float PrintWholeAlignment( const Sequence &query, const Sequence &target, const Cigar &cigar_ ) {
  std::string q;
  std::string t;
  std::string a;

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
          a += '|';
          q += query[ qcount++ ];
          t += target[ tcount++ ];
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
  std::cout << "QRY " << std::string( 11, ' ' ) << query.identifier << std::endl;
  std::cout << std::endl;
  std::cout << std::setw( 15 ) << queryStart + 1 << " " << q << " " << qcount << std::endl;
  std::cout << std::string( 16, ' ' ) << a << std::endl;
  std::cout << std::setw( 15 ) << targetStart + 1 << " " << t << " " << tcount << std::endl;
  std::cout << std::endl;
  std::cout << "REF " << std::string( 11, ' ' ) << target.identifier << std::endl;

  std::cout << std::endl;
  float identity = float( numMatches ) / float( numCols );
  std::cout <<  numCols << " cols, " << numMatches << " ids (" << ( 100.0f * identity ) << "%)" << std::endl;

  std::cout << std::endl;
  std::cout << std::string( 50, '=') << std::endl;

  return identity;
}


class Database {
  using SequenceRef = std::shared_ptr< Sequence >;
  using SequenceInfo = std::pair< size_t, SequenceRef >;
  using SequenceMappingDatabase = std::unordered_map< size_t, std::deque< SequenceInfo > >;

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
  Database( size_t wordSize )
    : mWordSize( wordSize )
  {
  }

  void AddSequence( const Sequence &seq ) {
    // Save
    SequenceRef ref = std::make_shared< Sequence >( seq );
    mSequences.push_back( ref );

    // Kmers for Indexing
    HashWords kmers( *ref, mWordSize );
    kmers.ForEach( [&]( size_t pos, size_t word ) {
      /* std::cout << pos << " " << word << std::endl; */
      this->mWordDB[ word ].push_back( SequenceInfo( pos, ref ) );
    });
  }

  SequenceList Query( const Sequence &query, float minIdentity, int maxHits = 1 ) {
    const size_t defaultMinHSPLength = 16;
    const size_t maxHSPJoinDistance = 16;

    size_t minHSPLength = std::min( defaultMinHSPLength, query.Length() / 2 );

    std::unordered_map< SequenceRef, HitTracker > hits;

    // Go through each kmer, find hits
    HashWords kmers( query, mWordSize );
    kmers.ForEach( [&]( size_t pos, size_t word ) {
      for( auto &seqInfo : mWordDB[ word ] ) {
        size_t candidatePos = seqInfo.first;
        SequenceRef candidateRef = seqInfo.second;
        hits[ candidateRef ].AddHit( pos, candidatePos, mWordSize );
      }
    });

    // Order candidates by number of shared words (hits)
    using Highscore = std::multimap< HitTracker, SequenceRef >;
    Highscore highscore;
    for( auto &h : hits ) {
      highscore.insert( std::make_pair( h.second, h.first ) );
    }

    // For each candidate:
    // - Get HSPs,
    // - Check for good HSP (>= similarity threshold)
    // - Join HSP together
    // - Align
    // - Check similarity
    /* std::cout << "=====NEWQUERY=====" << std::endl; */
    /* std::cout << "MinHSP Length " << minHSPLength << std::endl; */
    size_t numHits = 0;
    for( auto it = highscore.rbegin(); it != highscore.rend(); ++it ) {
      const HitTracker &hitTracker = it->first;
      const Sequence &candidateSeq = *( it->second );

      // some yolo preliminary optimization
      /* if( candidate.second.Score() < 0.5 * query.Length() ) { */
      /*   continue; */
      /* } */

      // Find all HSP
      // Sort by length
      // Try to find best chain
      // Fill space between with banded align

      std::set< HSP > hsps;

      for( auto &seed : hitTracker.Seeds() ) {
        Cigar leftCigar;
        size_t leftQuery, leftCandidate;

        int leftScore = mExtendAlign.Extend( query, candidateSeq,
            &leftQuery, &leftCandidate,
            &leftCigar,
            AlignmentDirection::backwards,
            seed.s1, seed.s2 );

        Cigar rightCigar;
        size_t rightQuery, rightCandidate;
        int rightScore = mExtendAlign.Extend( query, candidateSeq,
            &rightQuery, &rightCandidate,
            &rightCigar,
            AlignmentDirection::forwards,
            seed.s1 + seed.length, seed.s2 + seed.length );

        HSP hsp( leftQuery, rightQuery, leftCandidate, rightCandidate );
        if( hsp.Length() >= minHSPLength ) {
          // Construct hsp cigar
          // Middles matches 100% (exact seeds)
          Cigar middleCigar;
          middleCigar.push_front( CigarEntry( seed.length, CigarOp::MATCH ) );
          hsp.cigar = leftCigar + middleCigar + rightCigar ;

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
              next.a1 - 1, next.b1 - 1);
          alignment += cigar;
        }

        // Align last HSP's end to whole sequences end
        auto &last = *chain.crbegin();
        alignment += last.cigar;
        mBandedAlign.Align( query, candidateSeq, &cigar, AlignmentDirection::forwards, last.a2 + 1, last.b2 + 1);
        alignment += cigar;

        float identity = CalculateIdentity( alignment );
        if( identity >= minIdentity ) {
          PrintWholeAlignment( query, candidateSeq, alignment );

          numHits++;
          if( numHits >= maxHits )
            break;
        }
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
