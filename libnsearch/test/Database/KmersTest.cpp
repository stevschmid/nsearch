#include <catch.hpp>

#include <nsearch/Database/Kmers.h>

#include <sstream>

Kmer Kmerify( const std::string &a ) {
  Kmer kmer = 0;
  size_t counter = 0;
  for( auto &ch : a ) {
    char val;
    switch( ch ) {
      case 'A': val = 0b00; break;
      case 'C': val = 0b01; break;
      case 'G': val = 0b11; break;
      default: val = 0b10; /// T, U
    }
    kmer |= val << counter;
    counter += 2;
  }
  return kmer;
}

TEST_CASE( "Kmers" )  {
  std::vector< Kmer > out;

  SECTION( "Default" ) {
    Kmers k( "TAGAGAGAG", 4 );
    k.ForEach( [&]( Kmer kmer, size_t ) {
      out.push_back( kmer );
    });
  }

  SECTION( "Kmer length exceeds sequence" ) {
    Kmers k( "ATG", 4 );
    k.ForEach( [&]( Kmer kmer, size_t ) {
      out.push_back( kmer );
    });

    REQUIRE( out.size() == 1 );
    REQUIRE( out.front() == Kmerify( "ATG" ) );
  }

  SECTION( "Ambiguous Nucleotides") {
    Kmers k( "ATNCGTAT", 3 );
    k.ForEach( [&]( Kmer kmer, size_t ) {
      out.push_back( kmer );
    });

    REQUIRE( out.size() == 3 );
    REQUIRE( out.front() == Kmerify( "CGT" ) );

    Kmers k2( "ATGNTTA", 3 );
    out.clear();
    k2.ForEach( [&]( Kmer kmer, size_t ) {
      out.push_back( kmer );
    });

    REQUIRE( out.size() == 2 );
    REQUIRE( out.front() == Kmerify( "ATG" ) );
    REQUIRE( out.back() == Kmerify( "TTA" ) );
  }
}
