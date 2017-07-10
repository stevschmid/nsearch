#pragma once

#include "Sequence.h"
#include "Utils.h"

class Kmers {
public:
  Kmers( const Sequence &sequence, size_t length )
    : mSequence( sequence ), mLength( length )
  {
  }

  size_t Length() const {
    return mLength;
  }

  void ForEach( const std::function< void( const Sequence&, size_t pos ) > &block ) const {
    ForEveryKmerInSequence( [&]( size_t pos, Sequence &kmer ) {
      block( kmer, pos );
    });
  }

private:
  size_t mLength;
  Sequence mSequence;

  void ForEveryKmerInSequence( const std::function< void( size_t, Sequence& ) > &block ) const {
    for( int i = 0; i <= mSequence.Length() - mLength; i++ ) {
      Sequence kmer = mSequence.Subsequence( i, mLength );
      ForEveryKmerVariation( block, kmer, i );
    }
  }

  void ForEveryKmerVariation( const std::function< void( size_t, Sequence& ) > &block,
      Sequence &kmer, size_t global, size_t start = 0 ) const
  {
    static const char BASES[] = "ATCGU";

    // Find first ambiguous nucleotide
    size_t pos = kmer.sequence.find_first_not_of( BASES, start );
    if( pos != std::string::npos ) {
      bool skip = false;

      if( !skip ) {
        // Iterate through all possible bases
        char nuc = kmer[ pos ];
        for( int i = 0; i < sizeof( BASES ) - 1; i++ ) {
          if( DoNucleotidesMatch( nuc, BASES[ i ] ) )  {
            kmer.sequence[ pos ] = BASES[ i ];
            ForEveryKmerVariation( block, kmer, global, pos + 1 );
          }
        }
      }
    } else {
      block( global, kmer );
    }
  }
};
