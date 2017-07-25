#include "nsearch/Database/GlobalSearch.h"

#include "nsearch/Sequence.h"
#include "nsearch/Database.h"

#include "nsearch/Alignment/Common.h"
#include "nsearch/Alignment/Ranges.h"
#include "nsearch/Alignment/HitTracker.h"

GlobalSearch::GlobalSearch( const Database &db, float minIdentity, int maxHits, int maxRejects )
  : Search( db ), mDB( db ), mMinIdentity( minIdentity ), mMaxHits( maxHits ), mMaxRejects( maxRejects )
{
}

Search::ResultList GlobalSearch::Query( const Sequence &query )
{
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

  Highscore highscore( mMaxHits + mMaxRejects );

  Kmers kmers( query, mDB.mWordSize );
  std::vector< bool > uniqueCheck( mDB.mMaxUniqueWords );

  kmers.ForEach( [&]( Kmer word, size_t pos ) {
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

    kmers.ForEach( [&]( Kmer word, size_t pos ) {
        const Database::WordEntry *ptr = &mDB.mFirstEntries[ mDB.mIndexByWord[ word ] ];
        for( uint32_t i = 0; i < mDB.mNumEntriesByWord[ word ]; i++, ptr++ ) {
        if( ptr->sequence != seqIdx )
        continue;

        while( ptr ) {
        hitTracker.AddHit( pos, ptr->pos, kmers.Length() );
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
      if( identity >= mMinIdentity ) {
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
      if( numHits >= mMaxHits )
        break;
    } else {
      numRejects++;
      if( numRejects >= mMaxRejects )
        break;
    }
  }

  return Search::ResultList();
}
