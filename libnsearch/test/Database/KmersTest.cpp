#include <catch.hpp>

#include <nsearch/Database/Kmers.h>

#include <vector>
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
  Sequence seq;
  std::vector< Kmer > out;

  SECTION( "Default" ) {
    seq = "TAGAGAGAG";
    Kmers k( seq, 4 );
    k.ForEach( [&]( Kmer kmer, size_t ) {
      out.push_back( kmer );
    });
  }

  SECTION( "Kmer length exceeds sequence" ) {
    seq = "ATG";
    Kmers k( seq, 4 );
    k.ForEach( [&]( Kmer kmer, size_t ) {
      out.push_back( kmer );
    });

    REQUIRE( out.size() == 1 );
    REQUIRE( out.front() == Kmerify( "ATG" ) );
  }

  SECTION( "Ambiguous Nucleotides") {
    seq = "ATNCGTAT";
    Kmers k( seq, 3 );
    k.ForEach( [&]( Kmer kmer, size_t ) {
      out.push_back( kmer );
    });

    REQUIRE( out.size() == 3 );
    REQUIRE( out.front() == Kmerify( "CGT" ) );

    seq = "ATGNTTA";
    Kmers k2( seq, 3 );
    out.clear();
    k2.ForEach( [&]( Kmer kmer, size_t ) {
      out.push_back( kmer );
    });

    REQUIRE( out.size() == 2 );
    REQUIRE( out.front() == Kmerify( "ATG" ) );
    REQUIRE( out.back() == Kmerify( "TTA" ) );
  }
}
