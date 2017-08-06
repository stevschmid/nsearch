#pragma once

#include "../Sequence.h"

struct Protein {
  typedef char CharType;
};

// Based on BLOSUM62
// Collapse AAs into 4 bits
template <>
struct BitMapPolicy< Protein > {
  static const size_t NumBits = 4;

  inline static int8_t BitMap( const char aa ) {
    static const char BitMapping[] = {
      0b00000, // 'A'
      0b10000, // 'B' ambiguous/invalid
      0b00011, // 'C'
      0b00100, // 'D'
      0b00100, // 'E'
      0b01111, // 'F'
      0b00101, // 'G'
      0b00110, // 'H'
      0b00111, // 'I'
      0b10000, // 'J' ambiguous/invalid
      0b01001, // 'K'
      0b01000, // 'L'
      0b01010, // 'M'
      0b00010, // 'N'
      0b10000, // 'O' ambiguous/invalid
      0b01011, // 'P'
      0b00100, // 'Q'
      0b00001, // 'R'
      0b01100, // 'S'
      0b01101, // 'T'
      0b10000, // 'U' ambiguous/invalid
      0b00111, // 'V'
      0b01110, // 'W'
      0b10010, // 'X' ambiguous/invalid
      0b01111, // 'Y'
      0b10001, // 'Z' ambiguous/invalid
    };

    if( aa & 0b10010 > 0 )
      return -1;

    return BitMapping[ aa - 'A' ];
  }
};

template <>
struct ComparePolicy< Protein > {
  inline static int8_t Score( const char aaA, const char aaB ) {
    static const int ScoreMatrixSize = 26; // 'A'...'Z', waste some space for faster lookup
    static const int8_t ScoreMatrix[ ScoreMatrixSize ][ ScoreMatrixSize ] = {
    };

    return ScoreMatrix[ aaA - 'A' ][ aaB - 'A' ];
  }

  inline static bool Match( const char aaA, const char aaB ) {
    return Score( aaA, aaB ) > 0;
  }
};
