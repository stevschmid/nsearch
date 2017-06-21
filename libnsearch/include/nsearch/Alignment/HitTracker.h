#pragma once

#include "Seed.h"
#include "Ranges.h"

#include <map>

class HitTracker
{
  typedef std::pair< size_t, size_t > Pos;

public:
  void AddHit( size_t start1, size_t start2, size_t length ) {
    size_t diagStart1 = 0;
    size_t diagStart2 = start2 - start1;
    Ranges &ranges = mDiagonals[ Pos( diagStart1, diagStart2 ) ];
    ranges.Add( start1, start1 + length );
  }

  SeedList Seeds() const {
    SeedList seeds;

    for( auto &it : mDiagonals ) {
      auto &key = it.first;
      auto &ranges = it.second;

      for( auto &range : ranges.Get() ) {
        // (6,1) -> (8,3)
        size_t length = range.second - range.first;
        size_t start1 = range.first;
        size_t start2 = key.second + start1;
        seeds.push_back( Seed( start1, start2, length ) );
      }
    }

    return seeds;
  };

private:
  std::map< Pos, Ranges > mDiagonals;
};
