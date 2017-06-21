#pragma once

#include <map>
#include <cassert>

class RangeMerger {
public:
  size_t NumRanges() const {
    return mRanges.size();
  }

  size_t Size() const {
    size_t size = 0;
    for( auto &range : mRanges ) {
      size += ( range.second - range.first ) + 1;
    }
    return size;
  }

  // [rangeMin, rangeMax]
  void Add( size_t rangeMin, size_t rangeMax ) {
    assert( rangeMin <= rangeMax );

    std::map< size_t, size_t >::iterator insertIt, afterIt = mRanges.upper_bound( rangeMin );

    if( afterIt == mRanges.begin() || std::prev( afterIt )->second < rangeMin ) {
      insertIt = mRanges.insert( afterIt, std::pair< size_t, size_t >( rangeMin, rangeMax ) );
    } else {
      insertIt = std::prev( afterIt );
      if( insertIt->second < rangeMax ) {
        // Extend overlap
        insertIt->second = rangeMax;
      }
    }

    // Merge following ranges if covered fully by this one
    while( afterIt != mRanges.end() && rangeMax >= afterIt->first ) {
      insertIt->second = std::max( afterIt->second, insertIt->second );
      afterIt = mRanges.erase( afterIt );
    }
  }

  bool operator<( const RangeMerger &other ) const {
    return Size() < other.Size();
  }

  const std::map< size_t, size_t >& Ranges() const {
    return mRanges;
  }

private:
  std::map< size_t, size_t > mRanges; // Non overlapping ranges
};
