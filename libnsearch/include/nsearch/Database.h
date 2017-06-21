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

class Database {
  typedef std::pair< size_t, const Sequence* > SequenceInfo;

public:
  Database( int wordSize )
    : mWordSize( wordSize )
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

  SequenceList Query( const Sequence &query, int maxHits = 10 ) const {
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

    for( auto &c : candidates ) {
      OptimalChainFinder optimalChain( c.second.Seeds() );

      /* std::cout << "=======" << std::endl; */
      /* std::cout << (*c.first).identifier << std::endl; */
      /* std::cout << query.identifier << std::endl; */
      /* std::cout << (*c.first).Length() << std::endl; */
      /* for( auto &s : c.second.Seeds() ) { */
      /*   std::cout << s.s1 << " " << s.s2 << " " << s.length << std::endl; */
      /* } */
      /* std::cout << "Number of seeds" << c.second.Seeds().size() << std::endl; */
      /* std::cout << " Optimal chain size " << optimalChain.OptimalChain().size() << std::endl; */
    }

    // Sort candidates by coverage
    /* std::set< std::pair< Coverage, const Sequence* > > sortedCandidates; */
    /* for( auto &c : candidates ) */
    /*   sortedCandidates.insert( std::pair< Coverage, const Sequence* >( c.second, c.first ) ); */

    /* // Now score the candidates */
    SequenceList list;
    /* for( auto it = sortedCandidates.rbegin(); it != sortedCandidates.rend(); ++it ) { */
    /*   /1* std::cout << "Candidate " << (*it).first.CoveredFraction() << std::endl; *1/ */
    /*   list.push_back( *(*it).second ); */

    /*   if( list.size() >= maxHits ) */
    /*     break; */
    /* } */
    /* /1* std::cout << "====== " << std::endl; *1/ */

    return list;
  }

  size_t NumWords() const {
    return mWords.size();
  }

private:
  int mWordSize;
  SequenceList mSequences;
  std::unordered_map< Sequence, std::deque< SequenceInfo > > mWords;
};
