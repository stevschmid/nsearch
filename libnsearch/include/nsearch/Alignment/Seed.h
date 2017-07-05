#pragma once

#include <deque>

class Seed {
public:
  size_t s1, s2, length;

  Seed( size_t s1, size_t s2, size_t length )
    : s1( s1 ), s2( s2 ), length( length )
  {
  }

  bool operator==( const Seed &other ) const {
    return s1 == other.s1 && s2 == other.s2 && length == other.length;
  }
};

using SeedList = std::deque< Seed >;
