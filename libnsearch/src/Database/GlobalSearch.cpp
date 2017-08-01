#include "nsearch/Database/GlobalSearch.h"

#include "nsearch/Database.h"
#include "nsearch/Database/HSP.h"
#include "nsearch/Sequence.h"

#include "nsearch/Alignment/Common.h"

#include <set>

GlobalSearch::GlobalSearch( const Database& db, const float minIdentity,
                            const int maxHits, const int maxRejects )
    : BaseSearch( db ), mDB( db ), mMinIdentity( minIdentity ),
      mMaxHits( maxHits ), mMaxRejects( maxRejects ) {}

GlobalSearch::HitList GlobalSearch::Query( const Sequence& query ) {
  const size_t defaultMinHSPLength = 16;
  const size_t maxHSPJoinDistance  = 16;

  /* std::cout << "===> SEARCH " << query.identifier << std::endl; */

  size_t minHSPLength = std::min( defaultMinHSPLength, query.Length() / 2 );

  // Go through each kmer, find hits
  if( mHits.size() < mDB.Size() ) {
    mHits.resize( mDB.Size() );
  }

  // Fast counter reset
  memset( mHits.data(), 0, sizeof( Counter ) * mHits.capacity() );

  Highscore highscore( mMaxHits + mMaxRejects );

  auto hitsData = mHits.data();

  std::vector< Kmer > kmers;
  std::vector< bool > uniqueCheck( mDB.mMaxUniqueKmers, false );
  Kmers( query, mDB.mKmerLength ).ForEach( [&]( const Kmer kmer, const size_t pos ) {
    kmers.push_back( kmer );

    if( uniqueCheck[ kmer ] )
      return;

    uniqueCheck[ kmer ] = true;

    auto offset = mDB.mSequenceIdsOffsetByKmer[ kmer ];
    auto seqIds = &mDB.mSequenceIds[ offset ];
    for( auto i = 0; i < mDB.mSequenceIdsCountByKmer[ kmer ]; i++ ) {
      auto    seqId   = seqIds[ i ];
      Counter counter = ++hitsData[ seqId ];
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

  HitList hits;

  for( auto it = highscores.cbegin(); it != highscores.cend(); ++it ) {
    const size_t seqId = it->id;
    assert( seqId < mDB.mSequences.size() );
    const Sequence& candidateSeq = mDB.mSequences[ seqId ];

    std::deque< HSP > sps;

    for( size_t pos = 0; pos < kmers.size(); pos++ ) {
      auto kmers2 = &mDB.mKmers[ mDB.mKmerOffsetBySequenceId[ seqId ] ];
      for( size_t pos2 = 0; pos2 < mDB.mKmerCountBySequenceId[ seqId ];
           pos2++ ) {
        if( kmers2[ pos2 ] != kmers[ pos ] )
          continue;

        if( pos == 0 || pos2 == 0 ||
            ( kmers[ pos - 1 ] != kmers2[ pos2 - 1 ] ) ) {
          size_t length = mDB.mKmerLength;

          size_t cur  = pos + 1;
          size_t cur2 = pos2 + 1;
          while( cur < kmers.size() &&
                 cur2 < mDB.mKmerCountBySequenceId[ seqId ] &&
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
          bool match = DoNucleotidesMatch( query[ a ], candidateSeq[ b ] );
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
        hits.push_back( GlobalSearch::Hit{ candidateSeq, alignment } );
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
