#pragma once
#include "../Sequence.h"
#include "../Utils.h"

#define BIT_INDEX(pos) ( (pos) * 2 ) % ( sizeof( size_t ) * 8 )
#define BASE_VALUE(base) \
    ( (base) == 'A' ? 0b00 : \
      ( (base) == 'C' ? 0b01 : \
        ( (base) == 'T' || (base) == 'U' ) ? 0b10 : \
          ( (base) == 'G' ? 0b11 : (-1) ) \
      ) \
    )

using Kmer = uint32_t;

class Kmers {
public:
  using Callback = const std::function< void( Kmer, size_t ) >;

  Kmers( const Sequence &ref, size_t length )
    : mRef( ref )
  {
    mLength = std::min( { length, mRef.Length(), sizeof( Kmer ) * 4 } );
  }

  void ForEach( const Callback &block ) const {
    const char *ptr = mRef.sequence.data();

    // First kmer
    size_t lastAmbigIndex = ( size_t )-1;
    Kmer kmer = 0;
    for( size_t k = 0; k < mLength; k++ ) {
      int8_t val = BASE_VALUE( *ptr );
      if( val < 0 ) {
        lastAmbigIndex = k;
      } else {
        kmer |= ( val << BIT_INDEX( k ) );
      }
      ptr++;
    }

    if( lastAmbigIndex == ( size_t )-1 ) {
      block( kmer, 0 );
    } else {
      block( 0, 0 );
    }

    // For each consecutive kmer, shift window by one
    size_t maxFrame = mRef.Length() - mLength;
    for( size_t frame = 1; frame <= maxFrame; frame++, ptr++ ) {
      kmer >>= 2;
      int8_t val = BASE_VALUE( *ptr );
      if( val < 0 ) {
        lastAmbigIndex = frame + mLength - 1;
      } else {
        kmer |= ( val << BIT_INDEX( mLength - 1 ) );
      }

      if( lastAmbigIndex == ( size_t )-1 || frame > lastAmbigIndex ) {
        block( kmer, frame );
      } else {
        block( 0, frame );
      }
    }
  }

  size_t Count() const {
    return mRef.Length() - mLength + 1;
  }

private:
  size_t mLength;
  const Sequence &mRef;
};
