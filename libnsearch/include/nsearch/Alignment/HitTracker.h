#pragma once

#include "SegmentPair.h"
#include "Ranges.h"

#include <unordered_map>

class HitTracker
{
  using Diagonal = size_t;

public:
  void AddHit( size_t start1, size_t start2, size_t length ) {
    Diagonal diag = start2 - start1;
    Ranges &ranges = mDiagonals[ diag ];
    ranges.Add( start1, start1 + length );
  }

  SegmentPairList List() const {
    SegmentPairList segmentPairs;

    for( auto &it : mDiagonals ) {
      auto &diag = it.first;
      auto &ranges = it.second;

      for( auto &range : ranges.Get() ) {
        // (6,1) -> (8,3)
        size_t length = range.second - range.first;
        size_t start1 = range.first;
        size_t start2 = diag + start1;

        segmentPairs.emplace_back( start1, start2, length );
      }
    }

    return segmentPairs;
  };

private:
  std::unordered_map< Diagonal, Ranges > mDiagonals;
};
