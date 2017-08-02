#pragma once

#include "nsearch/Alignment/Cigar.h"
#include "nsearch/Sequence.h"
#include "nsearch/Database.h"

#include <deque>
#include <vector>

template< typename Alphabet >
struct Hit {
  Sequence< Alphabet > target;
  Cigar                alignment;
};

template< typename Alphabet >
using HitList = std::deque< Hit< Alphabet > >;

template< typename Alphabet >
class BaseSearch {
public:
  BaseSearch( const Database< Alphabet >& db ) : mDB( db ) {}
  virtual HitList< Alphabet > Query( const Sequence< Alphabet > & query ) = 0;

protected:
  const Database< Alphabet >& mDB;
};
