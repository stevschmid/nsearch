#pragma once

#include "../Sequence.h"

struct Protein {
  typedef char CharType;
};

template <>
struct BitMapPolicy< Protein > {
  static const size_t NumBits = 5; // 2^5 = 32 > 20

  inline static int8_t BitMap( const char aa ) {
    return aa - 'A';
  }
};
