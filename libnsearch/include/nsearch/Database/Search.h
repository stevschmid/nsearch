#pragma once

#include "nsearch/Alignment/Cigar.h"
#include "nsearch/Sequence.h"
#include "nsearch/Database.h"

#include "nsearch/Alphabet/DNA.h"

#include <deque>
#include <vector>

struct BaseSearchParams {
  int maxAccepts;
  int maxRejects;
  float minIdentity;
};

template < typename Alphabet >
struct SearchParams : public BaseSearchParams {
};

template<>
struct SearchParams< DNA > : public BaseSearchParams {
  DNA::Strand strand = DNA::Strand::Plus;
};

template < typename Alphabet >
struct Hit {
  Sequence< Alphabet > target;
  Cigar                alignment;
};

template<>
struct Hit< DNA > {
  Sequence< DNA > target;
  Cigar           alignment;
  DNA::Strand     strand;
};

template< typename Alphabet >
using HitList = std::deque< Hit< Alphabet > >;

template< typename Alphabet >
using SearchForHitsCallback =
  std::function< void( const Sequence< Alphabet >&, const Cigar& ) >;

template< typename Alphabet >
class Search {
public:
  Search( const Database< Alphabet >&     db,
          const SearchParams< Alphabet >& params )
      : mDB( db ), mParams( params ) {}

  inline HitList< Alphabet > Query( const Sequence< Alphabet >& query ) {
    HitList< Alphabet > hits;

    SearchForHits(
      query, [&]( const Sequence< Alphabet >& target, const Cigar& alignment ) {
        hits.push_back( { target, alignment } );
      } );

    return hits;
  }

protected:
  virtual void
  SearchForHits( const Sequence< Alphabet >&              query,
                 const SearchForHitsCallback< Alphabet >& callback ) = 0;

  const Database< Alphabet >&     mDB;
  const SearchParams< Alphabet >& mParams;
};

/*
 * For DNA, allow strand specification
 */
template <>
inline HitList< DNA > Search< DNA >::Query( const Sequence< DNA >& query ) {
  /* HitList< DNA > hits; */

  /* /1*   auto strand = mParams.strand; *1/ */

  /* /1*   if( strand == DNA::Strand::Plus || strand == DNA::Strand::Both ) { *1/ */
  /* /1*     /2* std::cout << "Search plus" << std::endl; *2/ *1/ */
  /* /1*     /2* auto ret = SearchForHits( query ); *2/ *1/ */
  /* /1*     SearchForHits( query, *1/ */
  /* /1*                    []( const Sequence< DNA >& target, const Cigar& */
  /*  * alignment ) { *1/ */
  /* /1*                      hits.push_back( *1/ */
  /* /1*                    } ); *1/ */
  /* /1*   } *1/ */

  /* /1* if( strand == DNA::Strand::Minus || strand == DNA::Strand::Both ) { *1/ */
  /* /1*   std::cout << "Search minus" << std::endl; *1/ */
  /* /1*   auto ret = SearchForHits( query.Reverse().Complement() ); *1/ */
  /* /1*   for( *1/ */
  /* /1* } *1/ */

  return HitList< DNA >();
}
