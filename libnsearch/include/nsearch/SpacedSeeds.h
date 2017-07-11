#pragma once

#include "Sequence.h"
#include "Utils.h"

#include <cstring>

#define BIT_INDEX(pos) ( (pos) * 2 ) % ( sizeof( size_t ) * 8 )
#define BASE_VALUE(base) \
    ( (base) == 'A' ? 0b00 : \
      ( (base) == 'C' ? 0b01 : \
        ( (base) == 'T' || (base) == 'U' ) ? 0b10 : \
          ( (base) == 'G' ? 0b11 : (-1) ) \
      ) \
    )

class SpacedSeeds {
public:
  using Callback = const std::function< void( size_t, size_t ) >;

public:
  SpacedSeeds( const Sequence &ref, size_t wordSize )
    : mRef( ref )
  {
    mWordSize = std::min( wordSize, mRef.Length() );
  }

  void ForEach( const Callback &block ) const {
    const char *ptr = mRef.sequence.data();

    size_t j;
    size_t cols = mRef.Length() - mWordSize + 1;
    for( size_t i = 0; i < cols; i++ ) {

      size_t word = 0;
      size_t counter = 0;
      for( j = 0; j < mWordSize; j++ ) {
        if( j % 3 == 2 )
          continue;

        int8_t val = BASE_VALUE( ptr[ j ] );
        if( val < 0 )
          break; // ambiguous nuc

        word |= ( val << BIT_INDEX( counter ) );
        counter++;
      }

      // check if word is unambiguous
      if( j == mWordSize )
        block( i, word );

      ptr++;
    }
  }

private:
  size_t mWordSize;
  const Sequence &mRef;
  std::string mPattern;
};
