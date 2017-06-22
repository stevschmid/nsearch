#pragma once

#include <unordered_map>
#include <set>

#include "Sequence.h"
#include "Utils.h"
#include "Aligner.h"
#include "Kmer.h"

#include "Alignment/Seed.h"
#include "Alignment/Ranges.h"
#include "Alignment/HitTracker.h"
#include "Alignment/OptimalChainFinder.h"
#include "Alignment/DP.h"

class Database {
  typedef std::pair< size_t, const Sequence* > SequenceInfo;

public:
  Database( int wordSize )
    : mWordSize( wordSize ), mDP( 16 )
  {
    // Sequences have to fit in the calculated hash
    assert( wordSize * 2 <= sizeof( size_t ) * 8 );
  }

  void AddSequence( const Sequence &seq ) {
    // Save
    mSequences.push_back( seq );
    const Sequence *ref = &mSequences.back();

    // Map words
    Kmers kmers( *ref, mWordSize );
    kmers.ForEach( [&]( const Sequence &kmer, size_t pos ) {
      this->mWords[ kmer ].push_back( SequenceInfo( pos, ref ) );
    });
  }

  SequenceList Query( const Sequence &query, int maxHits = 10 ) {
    std::unordered_map< const Sequence*, HitTracker > candidates;

    // Go through each kmer, find candidates
    Kmers kmers( query, mWordSize );
    kmers.ForEach( [&]( const Sequence &kmer, size_t pos ) {
      for( auto &seqInfo : mWords.at( kmer ) ) {
        size_t candidatePos = seqInfo.first;
        const Sequence *candidateSeq = seqInfo.second;

        candidates[ candidateSeq ].AddHit( pos, candidatePos, mWordSize );
      }
    });

    // Compute the optimal chain for each candidate
    std::multimap< OptimalChainFinder, const Sequence* > highscore;

    // Sort candidates based on the score of the optimal chain
    for( auto &c : candidates ) {
      const Sequence *seq = c.first;
      const HitTracker &hitTracker = c.second;

      // Compute optimal chain
      OptimalChainFinder ocf( hitTracker.Seeds() );
      highscore.insert( std::pair< OptimalChainFinder, const Sequence* >( ocf, seq ) );
    }

    Alignment aln;
    SequenceList list;
    std::cout << "===" << std::endl;
    std::cout << "QRY " << query.identifier << std::endl;

    for( auto it = highscore.rbegin(); it != highscore.rend(); ++it ) {
      const OptimalChainFinder &ocf = it->first;
      const Sequence &reference = *it->second;

      int score = mDP.AlignAlongChain( query, reference, ocf.OptimalChain() );
      list.push_back( *(*it).second );
      std::cout << reference.identifier << std::endl;
      std::cout << " Chain Score " << ocf.Score()
        << std::endl << " Align Score: " << score
        << std::endl << " Ref Length " << reference.Length() << std::endl;
      /* mDP.DebugPrint(); */

      if( list.size() >= maxHits )
        break;
    }
    std::cout << "=====" << std::endl;

    return list;
  }

  size_t NumWords() const {
    return mWords.size();
  }

private:
  GuidedBandedGlobalAlign mDP;
  int mWordSize;
  SequenceList mSequences;
  std::unordered_map< Sequence, std::deque< SequenceInfo > > mWords;
};
