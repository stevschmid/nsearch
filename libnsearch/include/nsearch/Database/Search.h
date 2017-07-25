#pragma once

#include "nsearch/Sequence.h"
#include "nsearch/Alignment/Cigar.h"

#include <deque>
#include <vector>

class Database;

class Search {
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

  Search( const Database &db ) : mDB( db ) { }
  virtual ResultList Query( const Sequence &query ) = 0;

private:
  const Database &mDB;
};

