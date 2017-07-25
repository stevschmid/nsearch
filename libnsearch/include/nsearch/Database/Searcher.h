#pragma once

#include "nsearch/Sequence.h"
#include "nsearch/Alignment/Cigar.h"
#include "nsearch/Alignment/ExtendAlign.h"
#include "nsearch/Alignment/BandedAlign.h"

#include <deque>
#include <vector>

class Database;

class DatabaseSearcher {
public:
  using Result = struct {
    Sequence query;
    Sequence target;

    Cigar alignment;

    size_t targetStart, targetLength;
    size_t queryStart, queryLength;

    size_t numCols, numMatches, numMismatches, numGaps;

    float identity;
  };
  using ResultList = std::deque< Result >;

  DatabaseSearcher( const Database &db );
  ResultList Query( const Sequence &query, float minIdentity, int maxHits = 1, int maxRejects = 8 );

private:
  const Database &mDB;

  std::vector< size_t > mHits;

  ExtendAlign mExtendAlign;
  BandedAlign mBandedAlign;
};
