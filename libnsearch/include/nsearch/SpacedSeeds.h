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
    static const std::string defaultPattern = "11101001010011011100";
    wordSize = std::min( wordSize, mRef.Length() );

    mPattern = "";
    size_t i = 0, numOnes = 0;
    while( numOnes < wordSize ) {
      char ch = defaultPattern[ i % defaultPattern.size() ];
      mPattern += ch;
      if( ch == '1' ) {
        numOnes++;
      }
      i++;
    }
    mMaxWordLength = i;
  }

  void ForEach( const Callback &block ) const {
    const char *ptr = mRef.sequence.data();

    size_t wordLength = std::min( mMaxWordLength, mRef.Length() );
    size_t cols = mRef.Length() - wordLength + 1;

    size_t j;
    for( int i = 0; i < cols; i++ ) {
      size_t word = 0;
      size_t counter = 0;
      for( j = 0; j < wordLength; j++ ) {
        if( mPattern[ j ] == '0' )
          continue;

        int8_t val = BASE_VALUE( ptr[ j ] );
        if( val < 0 )
          break; // ambiguous nuc

        word |= ( val << BIT_INDEX( counter ) );
        counter++;
      }

      // check if word is unambiguous
      if( j == wordLength )
        block( i, word );

      ptr++;
    }
  }

private:
  size_t mMaxWordLength;
  const Sequence &mRef;
  std::string mPattern;
};
