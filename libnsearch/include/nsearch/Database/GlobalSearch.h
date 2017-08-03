#pragma once

#include "BaseSearch.h"

#include "../Database.h"
#include "../Alignment/Common.h"
#include "../Alignment/BandedAlign.h"
#include "../Alignment/ExtendAlign.h"

#include <set>

using Counter = unsigned short;

template < typename Alphabet >
class GlobalSearch : public BaseSearch< Alphabet > {
public:
  GlobalSearch( const Database< Alphabet >& db, const float minIdentity,
                const int maxHits = 1, const int maxRejects = 8 );
  HitList< Alphabet > Query( const Sequence< Alphabet >& query );

private:
  using BaseSearch< Alphabet >::mDB;

  float mMinIdentity;
  int   mMaxHits;
  int   mMaxRejects;

  std::vector< Counter >  mHits;
  ExtendAlign< Alphabet > mExtendAlign;
  BandedAlign< Alphabet > mBandedAlign;
};

template < typename A >
GlobalSearch< A >::GlobalSearch( const Database< A >& db,
                                 const float minIdentity, const int maxHits,
                                 const int maxRejects )
    : BaseSearch< A >( db ), mMinIdentity( minIdentity ), mMaxHits( maxHits ),
      mMaxRejects( maxRejects ) {}

template < typename A >
HitList< A > GlobalSearch< A >::Query( const Sequence< A >& query ) {
  const size_t defaultMinHSPLength = 16;
  const size_t maxHSPJoinDistance  = 16;

  /* std::cout << "===> SEARCH " << query.identifier << std::endl; */

  size_t minHSPLength = std::min( defaultMinHSPLength, query.Length() / 2 );

  // Go through each kmer, find hits
  if( mHits.size() < mDB.NumSequences() ) {
    mHits.resize( mDB.NumSequences() );
  }

  // Fast counter reset
  memset( mHits.data(), 0, sizeof( Counter ) * mHits.capacity() );

  Highscore highscore( mMaxHits + mMaxRejects );

  auto hitsData = mHits.data();

  std::vector< Kmer > kmers;
  std::vector< bool > uniqueCheck( mDB.MaxUniqueKmers(), false );
  Kmers< A >( query, mDB.KmerLength() )
    .ForEach( [&]( const Kmer kmer, const size_t pos ) {
      kmers.push_back( kmer );

      if( kmer == AmbiguousKmer || uniqueCheck[ kmer ] )
        return;

      uniqueCheck[ kmer ] = true;

      size_t numSeqIds;
      const SequenceId* seqIds;

      if( !mDB.GetSequenceIdsIncludingKmer( kmer, &seqIds, &numSeqIds ) )
        return;

      for( size_t i = 0; i < numSeqIds; i++ ) {
        const auto& seqId   = seqIds[ i ];
        Counter     counter = ++hitsData[ seqId ];

        highscore.Set( seqId, counter );
      }
    } );

  // For each candidate:
  // - Get HSPs,
  // - Check for good HSP (>= similarity threshold)
  // - Join HSP together
  // - Align
  // - Check similarity
  int numHits    = 0;
  int numRejects = 0;

  auto highscores = highscore.EntriesFromTopToBottom();

  HitList< A > hits;

  for( auto it = highscores.cbegin(); it != highscores.cend(); ++it ) {
    const size_t    seqId        = it->id;
    const Sequence< A >& candidateSeq = mDB.GetSequenceById( seqId );

    std::deque< HSP > sps;

    for( size_t pos = 0; pos < kmers.size(); pos++ ) {
      const Kmer* kmers2;
      size_t      kmers2count;
      if( !mDB.GetKmersForSequenceId( seqId, &kmers2, &kmers2count ) )
        continue;

      for( size_t pos2 = 0; pos2 < kmers2count; pos2++ ) {
        if( kmers2[ pos2 ] != kmers[ pos ] )
          continue;

        if( pos == 0 || pos2 == 0 || kmers[ pos - 1 ] == AmbiguousKmer ||
            kmers2[ pos2 - 1 ] == AmbiguousKmer ||
            ( kmers[ pos - 1 ] != kmers2[ pos2 - 1 ] ) ) {
          size_t length = mDB.KmerLength();

          size_t cur  = pos + 1;
          size_t cur2 = pos2 + 1;
          while( cur < kmers.size() && cur2 < kmers2count &&
                 kmers[ cur ] != AmbiguousKmer &&
                 kmers2[ cur ] != AmbiguousKmer &&
                 kmers[ cur ] == kmers2[ cur2 ] ) {
            cur++;
            cur2++;
            length++;
          }

          sps.emplace_back( pos, cur - 1, pos2, cur2 - 1 );
        }
      }
    };

    // Find all HSP
    // Sort by length
    // Try to find best chain
    // Fill space between with banded align

    std::set< HSP > hsps;
    for( auto& sp : sps ) {
      size_t queryPos, candidatePos;

      size_t a1 = sp.a1, a2 = sp.a2, b1 = sp.b1, b2 = sp.b2;

      Cigar leftCigar;
      int   leftScore =
        mExtendAlign.Extend( query, candidateSeq, &queryPos, &candidatePos,
                             &leftCigar, AlignmentDirection::rev, a1, b1 );
      if( !leftCigar.empty() ) {
        a1 = queryPos;
        b1 = candidatePos;
      }

      Cigar  rightCigar;
      size_t rightQuery, rightCandidate;
      int    rightScore = mExtendAlign.Extend(
        query, candidateSeq, &queryPos, &candidatePos, &rightCigar,
        AlignmentDirection::fwd, a2 + 1, b2 + 1 );
      if( !rightCigar.empty() ) {
        a2 = queryPos;
        b2 = candidatePos;
      }

      HSP hsp( a1, a2, b1, b2 );
      if( hsp.Length() >= minHSPLength ) {
        // Construct hsp cigar (spaced seeds so we cannot assume full match)
        Cigar middleCigar;
        int   middleScore = 0;
        for( size_t a = sp.a1, b = sp.b1; a <= sp.a2 && b <= sp.b2; a++, b++ ) {
          bool match = A::Match( query[ a ], candidateSeq[ b ] );
          middleCigar.Add( match ? CigarOp::MATCH : CigarOp::MISMATCH );
          middleScore += match ? mExtendAlign.AP().matchScore
                               : mExtendAlign.AP().mismatchScore;
        }
        hsp.score = leftScore + middleScore + rightScore;
        hsp.cigar = leftCigar + middleCigar + rightCigar;

        // Save HSP
        hsps.insert( hsp );
      }
    }

    // Greedy join HSPs if close
    struct HSPChainOrdering {
      bool operator()( const HSP& left, const HSP& right ) const {
        return left.a1 < right.a1 && left.b1 < right.b1;
      }
    };

    std::set< HSP, HSPChainOrdering > chain;
    for( auto it = hsps.rbegin(); it != hsps.rend(); ++it ) {
      const HSP& hsp = *it;
      bool       hasNoOverlaps =
        std::none_of( chain.begin(), chain.end(), [&]( const HSP& existing ) {
          return hsp.IsOverlapping( existing );
        } );
      if( hasNoOverlaps ) {
        bool anyHSPJoinable =
          std::any_of( chain.begin(), chain.end(), [&]( const HSP& existing ) {
            return hsp.DistanceTo( existing ) <= maxHSPJoinDistance;
          } );

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
      auto& first = *chain.cbegin();
      mBandedAlign.Align( query, candidateSeq, &cigar, AlignmentDirection::rev,
                          first.a1, first.b1 );
      alignment += cigar;

      // Align in between the HSP's
      for( auto it1 = chain.cbegin(), it2 = ++chain.cbegin();
           it1 != chain.cend() && it2 != chain.cend(); ++it1, ++it2 ) {
        auto& current = *it1;
        auto& next    = *it2;

        alignment += current.cigar;
        mBandedAlign.Align( query, candidateSeq, &cigar,
                            AlignmentDirection::fwd, current.a2 + 1,
                            current.b2 + 1, next.a1, next.b1 );
        alignment += cigar;
      }

      // Align last HSP's end to whole sequences end
      auto& last = *chain.crbegin();
      alignment += last.cigar;
      mBandedAlign.Align( query, candidateSeq, &cigar, AlignmentDirection::fwd,
                          last.a2 + 1, last.b2 + 1 );
      alignment += cigar;

      float identity = alignment.Identity();
      if( identity >= mMinIdentity ) {
        accept = true;
        hits.push_back( { candidateSeq, alignment } );
      }
    }

    if( accept ) {
      numHits++;
      if( numHits >= mMaxHits )
        break;
    } else {
      numRejects++;
      if( numRejects >= mMaxRejects )
        break;
    }
  }

  return hits;
}
