#pragma once

#include "Sequence.h"
#include "Utils.h"

#include <cstring>
#include <tmmintrin.h>

#define SHL128(v, n) \
({ \
    __m128i v1, v2; \
 \
    if ((n) >= 64) \
    { \
        v1 = _mm_slli_si128(v, 8); \
        v1 = _mm_slli_epi64(v1, (n) - 64); \
    } \
    else \
    { \
        v1 = _mm_slli_epi64(v, n); \
        v2 = _mm_slli_si128(v, 8); \
        v2 = _mm_srli_epi64(v2, 64 - (n)); \
        v1 = _mm_or_si128(v1, v2); \
    } \
    v1; \
})

#define SHR128(v, n) \
({ \
    __m128i v1, v2; \
 \
    if ((n) >= 64) \
    { \
        v1 = _mm_srli_si128(v, 8); \
        v1 = _mm_srli_epi64(v1, (n) - 64); \
    } \
    else \
    { \
        v1 = _mm_srli_epi64(v, n); \
        v2 = _mm_srli_si128(v, 8); \
        v2 = _mm_slli_epi64(v2, 64 - (n)); \
        v1 = _mm_or_si128(v1, v2); \
    } \
    v1; \
 })

class SpacedSeedsSIMD {
public:
  using Callback = const std::function< void( size_t, size_t ) >;
  static const int PATTERN_SIZE = 16;

public:
  SpacedSeedsSIMD( const Sequence &ref, size_t numOnes )
    : mRef( ref )
  {
    static const char BASE_PATTERN[] = "111010010100110111";

    // 111010010100110111
    //   __m128i indexMask = _mm_setr_epi8(
    // 0, 1, 2, 4,
    //   7, 9, 12, 13,
    //   15, 0x80, 0x80, 0x80,
    //   0x80, 0x80, 0x80, 0x80
    //     );

    alignas(16) uint8_t indices[ PATTERN_SIZE ];
    memset( indices, 0x80, PATTERN_SIZE );

    int i, countOnes;
    for( i = 0, countOnes = 0;
        i < PATTERN_SIZE && countOnes < numOnes;
        i++ )
    {
      if( BASE_PATTERN[ i % sizeof( BASE_PATTERN ) ] == '1' ) {
        indices[ countOnes++ ] = i;
      }
    }
    mWordSize = i;
    mIndexMask = _mm_load_si128( (__m128i*)indices );
  }

  static bool ExtractSpacedSeed( const char *str, const char *end, __m128i maskIndex, uint32_t *seedOut ) {
    alignas(16) uint8_t buffer[16];
    memset( buffer, 'A', 16 );
    memcpy( buffer, str, end - str < 16 ? end - str : 16 );
    __m128i seq = _mm_load_si128( (__m128i*)buffer );

    // Convert U to T
    __m128i maskU = _mm_cmpeq_epi8( seq, _mm_set1_epi8( 'U' ) );
    seq = _mm_or_si128(
        _mm_andnot_si128( maskU, seq ),
        _mm_and_si128( maskU, _mm_set1_epi8( 'T' ) )
        );

    // Check for any ambiguous nucleotide
    __m128i maskA = _mm_cmpeq_epi8( seq, _mm_set1_epi8( 'A' ) );
    __m128i maskT = _mm_cmpeq_epi8( seq, _mm_set1_epi8( 'T' ) );
    __m128i maskC = _mm_cmpeq_epi8( seq, _mm_set1_epi8( 'C' ) );
    __m128i maskG = _mm_cmpeq_epi8( seq, _mm_set1_epi8( 'G' ) );

    __m128 maskATCG = _mm_or_si128(
        _mm_or_si128( maskA, maskT ),
        _mm_or_si128( maskC, maskG )
        );

    // Disregard ambiguous nucs which represent spaces in mask (ie 0)
    if( _mm_movemask_epi8( maskATCG ) != 0xFFFF ) {
      // Check if the ambiguous nuc is not space (shuffle is slow)
      __m128i check = _mm_or_si128(
          _mm_shuffle_epi8( maskATCG, maskIndex ),
          _mm_cmpeq_epi8( maskIndex, _mm_set1_epi8( 0x80 ) )
          );

      if( _mm_movemask_epi8( check ) != 0xFFFF ) {
        return false;
      }
    }

    // Convert to 2bits per base
    __m128i bin;
    bin = _mm_setzero_si128();
    bin = _mm_or_si128( bin, _mm_and_si128( maskC, _mm_set1_epi8( 0b01 ) ) );
    bin = _mm_or_si128( bin, _mm_and_si128( maskT, _mm_set1_epi8( 0b10 ) ) );
    bin = _mm_or_si128( bin, _mm_and_si128( maskG, _mm_set1_epi8( 0b11 ) ) );

    // "Remove" spaced nucs
    bin = _mm_shuffle_epi8( bin, maskIndex );

    // Bitshift to get 4 nucleotide into a single byte
    __m128i packed;
    packed = _mm_or_si128( bin, SHL128( _mm_bslli_si128( bin, 1 ), 2 ) );
    packed = _mm_or_si128( packed, SHL128( _mm_bslli_si128( packed, 1 ), 2 ) );
    packed = _mm_or_si128( packed, SHL128( _mm_bslli_si128( packed, 1 ), 2 ) );

    // Copy the final bytes
    __m128i final = _mm_shuffle_epi8( packed, _mm_set_epi8(
          0x80, 0x80, 0x80, 0x80,
          0x80, 0x80, 0x80, 0x80,
          0x80, 0x80, 0x80, 0x80,
          15, 11, 7, 3 ) );

    // Now convert
    *seedOut = _mm_cvtsi128_si32( final );
    return true;
  }

  void ForEach( const Callback &block ) const {
    const char *end = mRef.sequence.data() + mRef.sequence.length();
    const char *ptr = mRef.sequence.data();

    signed int cols = mRef.Length() - mWordSize + 1;
    for( signed int i = 0; i < cols; i++ ) {
      // check if word is unambiguous
      uint32_t word;
      if( ExtractSpacedSeed( ptr, end, mIndexMask, &word ) )
        block( i, word );

      ptr++;
    }
  }

private:
  const Sequence &mRef;
  size_t mWordSize;
  __m128i mIndexMask;
};
