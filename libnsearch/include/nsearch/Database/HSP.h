#pragma once

#include "../Alignment/Cigar.h"

// High-scoring segment pairss
// HSP: first and last character in sequence (i.e. seq[a1] - seq[a2])
class HSP {
public:
  size_t a1, a2;
  size_t b1, b2;
  int score;
  Cigar cigar;

  HSP( size_t a1, size_t a2, size_t b1, size_t b2, int score = 0 )
    : a1( a1 ), a2( a2 ), b1( b1 ), b2( b2 ), score( score )
  {
    assert( a2 >= a1 && b2 >= b1 );
  }

  size_t Length() const {
    return std::max( a2 - a1, b2 - b1 ) + 1;
  }

  int Score() const {
    return score;
  }

  bool IsOverlapping( const HSP &other ) const {
    return
      ( a1 <= other.a2 && other.a1 <= a2 ) // overlap in A direction
      ||
      ( b1 <= other.b2 && other.b1 <= b2 ); // overlap in B direction
  }

  size_t DistanceTo( const HSP &other ) const {
    size_t x = ( a1 > other.a2 ? a1 - other.a2 : other.a1 - a2 );
    size_t y = ( b1 > other.b2 ? b1 - other.b2 : other.b2 - b2 );

    return sqrt( x*x + y*y );
  }

  bool operator<( const HSP &other ) const {
    return Score() < other.Score();
  }
};
