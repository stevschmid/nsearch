#include <catch.hpp>

#include <nsearch/Database/Kmers.h>

#include <vector>

#include "../Support.h"

TEST_CASE( "Kmers" ) {
  Sequence            seq;
  std::vector< Kmer > out;

  SECTION( "Default" ) {
    seq = "TAGAGAGAG";
    Kmers k( seq, 4 );
    k.ForEach( [&]( Kmer kmer, size_t ) { out.push_back( kmer ); } );
  }

  SECTION( "Kmer length exceeds sequence" ) {
    seq = "ATG";
    Kmers k( seq, 4 );
    k.ForEach( [&]( Kmer kmer, size_t ) { out.push_back( kmer ); } );

    REQUIRE( out.size() == 1 );
    REQUIRE( out.front() == Kmerify( "ATG" ) );
  }

  SECTION( "Ambiguous Nucleotides" ) {
    seq = "ATNCGTAT";
    Kmers k( seq, 3 );
    k.ForEach( [&]( Kmer kmer, size_t ) { out.push_back( kmer ); } );

    REQUIRE( out.size() == 6 );
    REQUIRE( out.front() == InvalidKmer );

    seq = "ATGNTTA";
    Kmers k2( seq, 3 );
    out.clear();
    k2.ForEach( [&]( Kmer kmer, size_t ) { out.push_back( kmer ); } );

    REQUIRE( out.size() == 5 );
    REQUIRE( out[ 0 ] == Kmerify( "ATG" ) );
    REQUIRE( out[ 1 ] == InvalidKmer );
    REQUIRE( out[ 2 ] == InvalidKmer );
    REQUIRE( out[ 3 ] == InvalidKmer );
    REQUIRE( out[ 4 ] == Kmerify( "TTA" ) );
  }
}
