#pragma once

#include "nsearch/Alignment/Cigar.h"
#include "nsearch/Sequence.h"
#include "nsearch/Database.h"

#include <deque>
#include <vector>

template < typename Alphabet >
struct Hit {
  Sequence< Alphabet > target;
  Cigar                alignment;
};

struct BaseSearchParams {
  int maxAccepts;
  int maxRejects;
  float minIdentity;
};

template < typename Alphabet >
struct SearchParams : public BaseSearchParams {
};

template< typename Alphabet >
using HitList = std::deque< Hit< Alphabet > >;

template< typename Alphabet >
class Search {
public:
  Search( const Database< Alphabet >&     db,
          const SearchParams< Alphabet >& params )
      : mDB( db ), mParams( params ) {}

  HitList< Alphabet > Query( const Sequence< Alphabet >& query ) {
    return QuerySingleSequence( query );
  }

protected:
  virtual HitList< Alphabet > QuerySingleSequence( const Sequence< Alphabet > & query ) = 0;

  const Database< Alphabet >& mDB;
  const SearchParams< Alphabet >& mParams;
};
