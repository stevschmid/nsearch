#pragma once

#include <unordered_map>
#include <set>

#include "Sequence.h"
#include "Utils.h"
#include "Aligner.h"
#include "Kmer.h"

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
    std::unordered_map< const Sequence*, Coverage > candidates;

    // Go through each kmer, find candidates
    Kmers kmers( query, mWordSize );
    kmers.ForEach( [&]( const Sequence &kmer, size_t pos ) {
      for( auto &seqInfo : mWords.at( kmer ) ) {
        size_t candidatePos = seqInfo.first;
        const Sequence *candidateSeq = seqInfo.second;

        if( candidates.find( candidateSeq ) == candidates.end() ) {
          candidates[ candidateSeq ] = Coverage( query.Length() );
        }

        Coverage &coverage = candidates[ candidateSeq ];
        coverage.Add( pos, pos + mWordSize - 1 );
      }
    });

    // Sort candidates by coverage
    std::set< std::pair< Coverage, const Sequence* > > sortedCandidates;
    for( auto &c : candidates )
      sortedCandidates.insert( std::pair< Coverage, const Sequence* >( c.second, c.first ) );

    // Now score the candidates
    SequenceList list;
    for( auto it = sortedCandidates.rbegin(); it != sortedCandidates.rend(); ++it ) {
      /* std::cout << "Candidate " << (*it).first.CoveredFraction() << std::endl; */
      list.push_back( *(*it).second );

      if( list.size() >= maxHits )
        break;
    }
    /* std::cout << "====== " << std::endl; */

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
