#pragma once

#include <unordered_map>
#include <set>

#include "Sequence.h"
#include "Utils.h"
#include "Aligner.h"
#include "Kmer.h"

#include "Alignment/Alignment.h"
#include "Alignment/OptimalChainFinder.h"

class HitTracker
{
public:
  void AddHit( size_t start1, size_t start2, size_t length ) {
    size_t diagStart1 = 0;
    size_t diagStart2 = start2 - start1;
    Coverage &cov = mDiagonals[ std::pair< size_t, size_t >( diagStart1, diagStart2 ) ];
    cov.Add( start1, start1 + length );
  }

  SeedList Seeds() const {
    SeedList seeds;

    for( auto &it : mDiagonals ) {
      auto &key = it.first;
      auto &coverage = it.second;

      for( auto &range : coverage.Ranges() ) {
        // (6,1) -> (8,3)
        size_t length = range.second - range.first;
        size_t start1 = range.first;
        size_t start2 = key.second + start1;
        seeds.push_back( Seed( start1, start2, length ) );
      }
    }

    return seeds;
  };

private:
  std::map< std::pair< size_t, size_t >, Coverage > mDiagonals;
};

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

        /* if( candidates.find( candidateSeq ) == candidates.end() ) { */
        /*   candidates[ candidateSeq ] = Coverage( query.Length() ); */
        /* } */
        candidates[ candidateSeq ].AddHit( pos, candidatePos, mWordSize );
        /* Coverage &coverage = candidates[ candidateSeq ]; */
        /* coverage.Add( pos, pos + mWordSize - 1 ); */
      }
    });

    for( auto &c : candidates ) {
      if( c.second.Seeds().size() > 1 ) {
        if( (*c.first).identifier == "AF189721.1 Citrobacter amalonaticus beta-lactamase CTX-M-8 precursor (blaCTX-M-8) gene, complete cds"
            && query.identifier == "0-632-AY954516.1 Escherichia coli beta-lactamase CTX-M-39 gene, complete cds" )
          {
          OptimalChainFinder optimalChain( c.second.Seeds() );

          std::cout << "=======" << std::endl;
          std::cout << (*c.first).identifier << std::endl;
          std::cout << query.identifier << std::endl;
          std::cout << (*c.first).Length() << std::endl;
          for( auto &s : c.second.Seeds() ) {
            std::cout << s.s1 << " " << s.s2 << " " << s.length << std::endl;
          }
          std::cout << "Number of seeds" << c.second.Seeds().size() << std::endl;
          std::cout << " Optimal chain size " << optimalChain.OptimalChain().size() << std::endl;
        }
      }
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
