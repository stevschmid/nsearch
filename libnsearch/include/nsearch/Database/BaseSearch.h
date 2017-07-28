#pragma once

#include "nsearch/Sequence.h"
#include "nsearch/Alignment/Cigar.h"

#include <deque>
#include <vector>

class Database;

class BaseSearch {
public:
  using Hit = struct {
    Sequence target;
    Cigar alignment;
  };
  using HitList = std::deque< Hit >;

  BaseSearch( const Database &db ) : mDB( db ) { }
  virtual HitList Query( const Sequence &query ) = 0;

private:
  const Database &mDB;
};
