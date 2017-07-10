#pragma once

#include "Sequence.h"
#include "Utils.h"

#include <cstring>

#define BIT_INDEX(pos) ( (pos) * 2 ) % ( sizeof( size_t ) * 8 )
#define BASE_VALUE(base) \
    ( (base) == 'A' ? 0b00 : \
      ( (base) == 'C' ? 0b01 : \
        ( ( (base) == 'T' || (base) == 'U' ) ? 0b10 : 0b11 /* G */ ) \
      ) \
    )

class HashWords {
public:
  using Callback = const std::function< void( size_t, size_t ) >;

private:
  void ForEachVariation( const Callback &block, size_t pos, size_t word, size_t idx = 0 ) const {
    static const std::string BASES = "ACTGU";

    const char *ptr = mRef.sequence.data() + pos + idx;
    const char *cur = ptr;
    const char *end = mRef.sequence.data() + pos + mWordSize;

    while( cur < end ) {
      char ch = *cur;
      if( ch >= 97 && ch <= 122 ) // upcase
        ch &= ~0x20;
      if( ch != 'A' && ch != 'T' && ch != 'C' && ch != 'G' && ch != 'U' )
        break;
      cur++;
    }

    if( cur < end ) {
      size_t newIdx = idx + ( cur - ptr );
      for( auto &base : BASES ) {
        /* std::cout << *cur << std::endl; */
        if( DoNucleotidesMatch( *cur, base ) ) {
          word &= ~( 0b11 << BIT_INDEX( newIdx ) ); // clear
          word |= ( BASE_VALUE( base ) << BIT_INDEX( newIdx ) );
          ForEachVariation( block, pos, word, newIdx + 1 );
        }
      }
    } else {
      block( pos, word );
    }
  }

public:
  HashWords( const Sequence &ref, size_t wordSize )
    : mRef( ref )
  {
    mWordSize = std::min( wordSize, mRef.Length() );
  }

  void ForEach( const Callback &block ) const {
    const char *data = mRef.sequence.data();
    const char *ptr = data;

    // First word
    size_t word = 0;
    for( size_t k = 0; k < mWordSize; k++ ) {
      word |= ( BASE_VALUE( *ptr ) << BIT_INDEX( k ) );
      ptr++;
    }
    ForEachVariation( block, 0, word );

    // For each, shift window by one
    size_t cols = mRef.Length() - mWordSize;
    for( size_t i = 1; i <= cols; i++ ) {
      word >>= 2;
      word |= ( BASE_VALUE( *ptr ) << BIT_INDEX( mWordSize - 1 ) );
      ForEachVariation( block, i, word );
      ptr++;
    }
  }

private:
  size_t mWordSize;
  const Sequence &mRef;
};
