#pragma once
#include "../Sequence.h"
#include "../Utils.h"

using Kmer = uint32_t;

const Kmer AmbiguousKmer = ( Kmer )-1;

inline size_t BitIndexForPosition( const int pos ) {
  return ( pos * 2 ) % ( sizeof( Kmer ) * 8 );
}

inline int8_t MapNucleotide( const char base ) {
  switch( base ) {
    case 'A': return 0b00;
    case 'C': return 0b01;
    case 'U':
    case 'T': return 0b10;
    case 'G': return 0b11;
    default: return -1;
  }
}

class Kmers {
public:
  using Callback = const std::function< void( const Kmer, const size_t ) >;

  Kmers( const Sequence& ref, const size_t length ) : mRef( ref ) {
    mLength = std::min( { length, mRef.Length(), sizeof( Kmer ) * 4 } );
  }

  void ForEach( const Callback& block ) const {
    const char* ptr = mRef.sequence.data();

    // First kmer
    size_t lastAmbigIndex = ( size_t ) -1;
    Kmer   kmer           = 0;
    for( size_t k = 0; k < mLength; k++ ) {
      int8_t val = MapNucleotide( *ptr );
      if( val < 0 ) {
        lastAmbigIndex = k;
      } else {
        kmer |= ( val << BitIndexForPosition( k ) );
      }
      ptr++;
    }

    if( lastAmbigIndex == ( size_t ) -1 ) {
      block( kmer, 0 );
    } else {
      block( AmbiguousKmer, 0 );
    }

    // For each consecutive kmer, shift window by one
    size_t maxFrame = mRef.Length() - mLength;
    for( size_t frame = 1; frame <= maxFrame; frame++, ptr++ ) {
      kmer >>= 2;
      int8_t val = MapNucleotide( *ptr );
      if( val < 0 ) {
        lastAmbigIndex = frame + mLength - 1;
      } else {
        kmer |= ( val << BitIndexForPosition( mLength - 1 ) );
      }

      if( lastAmbigIndex == ( size_t ) -1 || frame > lastAmbigIndex ) {
        block( kmer, frame );
      } else {
        block( AmbiguousKmer, frame );
      }
    }
  }

  size_t Count() const {
    return mRef.Length() - mLength + 1;
  }

private:
  size_t          mLength;
  const Sequence& mRef;
};
