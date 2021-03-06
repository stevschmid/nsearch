#include <catch.hpp>

#include <nsearch/Database/Kmers.h>
#include <nsearch/Alphabet/DNA.h>

#include <vector>

#include "../Support.h"

TEST_CASE( "Kmers" ) {
  Sequence< DNA >     seq;
  std::vector< Kmer > out;

  SECTION( "Default" ) {
    seq = "TAGAGAGAG";
    Kmers< DNA > k( seq, 4 );
    k.ForEach( [&]( Kmer kmer, size_t ) { out.push_back( kmer ); } );
  }

  SECTION( "Kmer length exceeds sequence" ) {
    seq = "ATG";
    Kmers< DNA > k( seq, 4 );
    k.ForEach( [&]( Kmer kmer, size_t ) { out.push_back( kmer ); } );

    REQUIRE( out.size() == 1 );
    REQUIRE( out.front() == Kmerify( "ATG" ) );
  }

  SECTION( "Ambiguous Nucleotides" ) {
    seq = "ATNCGTAT";
    Kmers< DNA > k( seq, 3 );
    k.ForEach( [&]( Kmer kmer, size_t ) { out.push_back( kmer ); } );

    REQUIRE( out.size() == 6 );
    REQUIRE( out.front() == AmbiguousKmer );

    seq = "ATGNTTA";
    Kmers< DNA > k2( seq, 3 );
    out.clear();
    k2.ForEach( [&]( Kmer kmer, size_t ) { out.push_back( kmer ); } );

    REQUIRE( out.size() == 5 );
    REQUIRE( out[ 0 ] == Kmerify( "ATG" ) );
    REQUIRE( out[ 1 ] == AmbiguousKmer );
    REQUIRE( out[ 2 ] == AmbiguousKmer );
    REQUIRE( out[ 3 ] == AmbiguousKmer );
    REQUIRE( out[ 4 ] == Kmerify( "TTA" ) );
  }
}
