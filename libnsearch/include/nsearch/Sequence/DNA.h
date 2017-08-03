#pragma once

#include "../Sequence.h"
#include "../Utils.h"

struct DNA {
  typedef char CharType;
};

using RNA = DNA;

template <>
struct BitMapPolicy< DNA > {
  static const size_t NumBits = 2;

  inline static int8_t BitMap( const char base ) {
    switch( base ) {
      case 'A': return 0b00;
      case 'C': return 0b01;
      case 'U':
      case 'T': return 0b10;
      case 'G': return 0b11;
      default: return -1;
    }
  }
};

template <>
struct ComplementPolicy< DNA > {
  inline static char Complement( const char nuc ) {
    switch( nuc ) {
      case 'A': return 'T';
      case 'G': return 'C';
      case 'C': return 'G';
      case 'T': return 'A';
      case 'U': return 'A';

      case 'Y': return 'R';
      case 'R': return 'Y';
      case 'W': return 'W';
      case 'S': return 'S';
      case 'K': return 'M';
      case 'M': return 'K';

      case 'D': return 'H';
      case 'V': return 'B';
      case 'H': return 'D';
      case 'B': return 'V';
      case 'N': return 'N';
    }

    return nuc;
  }
};

template <>
struct MatchPolicy< DNA > {
  inline static bool Match( const char nucA, const char nucB ) {
    return DoNucleotidesMatch( nucA, nucB );
  }
};

