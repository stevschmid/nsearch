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
  using SequenceRef = std::shared_ptr< Sequence >;
  using SequenceInfo = std::pair< size_t, SequenceRef >;

public:
  Database( int wordSize )
    : mWordSize( wordSize ), mDP( 32 )
  {
    // Sequences have to fit in the calculated hash
    assert( wordSize * 2 <= sizeof( size_t ) * 8 );
  }

  void AddSequence( const Sequence &seq ) {
    // Save
    mSequences.push_back( std::make_shared< Sequence >( seq ) );
    SequenceRef ref = mSequences.back();

    // Map words
    Kmers kmers( *ref, mWordSize );
    kmers.ForEach( [&]( const Sequence &kmer, size_t pos ) {
      this->mWords[ kmer ].push_back( SequenceInfo( pos, ref ) );
    });
  }

  SequenceList Query( const Sequence &query, int maxHits = 10 ) {
    std::unordered_map< SequenceRef, HitTracker > candidates;

    // Go through each kmer, find candidates
    Kmers kmers( query, mWordSize );
    kmers.ForEach( [&]( const Sequence &kmer, size_t pos ) {
      auto mapIt = mWords.find( kmer );
      if( mapIt != mWords.end() ) {
        for( auto &seqInfo : mapIt->second ) {
          size_t candidatePos = seqInfo.first;
          SequenceRef candidateSeqRef = seqInfo.second;
          candidates[ candidateSeqRef ].AddHit( pos, candidatePos, mWordSize );
        }
      }
    });

    // Compute the optimal chain for each candidate
    std::multimap< OptimalChainFinder, SequenceRef > highscore;

    // Sort candidates based on the score of the optimal chain
    for( auto &c : candidates ) {
      SequenceRef seq = c.first;
      const HitTracker &hitTracker = c.second;

      // Compute optimal chain
      OptimalChainFinder ocf( hitTracker.Seeds() );
      highscore.insert( std::pair< OptimalChainFinder, SequenceRef >( ocf, seq ) );
    }

    Alignment aln;
    SequenceList list;

    for( auto it = highscore.rbegin(); it != highscore.rend(); ++it ) {
      const OptimalChainFinder &ocf = it->first;
      const Sequence &reference = *it->second;

      int score = mDP.AlignAlongChain( query, reference, ocf.OptimalChain(), &aln );
      list.push_back( *(*it).second );

      std::cout << "Query " << query.identifier << std::endl;
      std::cout << "Reference " << reference.identifier << std::endl;
      std::cout << aln << std::endl;
      std::cout << " Chain Score " << ocf.Score()
        << std::endl << " Align Score: " << score
        << std::endl << " Ref Length " << reference.Length() << std::endl;
      /* mDP.DebugPrint( true ); */

      if( list.size() >= maxHits )
        break;
    }

    return list;
  }

  size_t NumWords() const {
    return mWords.size();
  }

private:
  GuidedBandedGlobalAlign mDP;
  int mWordSize;
  std::deque< SequenceRef > mSequences;
  std::unordered_map< Sequence, std::deque< SequenceInfo > > mWords;
};
