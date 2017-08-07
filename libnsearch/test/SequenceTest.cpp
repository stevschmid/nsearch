#include <catch.hpp>

#include <nsearch/Sequence.h>
#include <nsearch/Alphabet/DNA.h>

TEST_CASE( "Sequence" ) { Sequence< DNA > seq( "Id", "ACCT", "JJ::" );

  REQUIRE( seq.Length() == 4 );
  REQUIRE( seq.identifier == "Id" );
  REQUIRE( seq.sequence == "ACCT" );
  REQUIRE( seq.quality == "JJ::" );

  SECTION( "subsequence" ) {
    Sequence< DNA > sub = seq.Subsequence( 2, 2 );
    REQUIRE( sub.sequence == "CT" );
    REQUIRE( sub.quality == "::" );

    REQUIRE( seq.Subsequence( 2 ).sequence == "CT" );
    REQUIRE( seq.Subsequence( 4 ).sequence == "" );

    Sequence< DNA > sub3 = Sequence< DNA >( "ATTG" ).Subsequence( 10, 10 );
    REQUIRE( sub3.sequence == "" );
  }

  SECTION( "operator +" ) {
    Sequence< DNA > m = seq + Sequence< DNA >( "Id2", "GCTT", "BBAJ" );

    REQUIRE( m.identifier == "Id" );
    REQUIRE( m.sequence == "ACCTGCTT" );
    REQUIRE( m.quality == "JJ::BBAJ" );
  }

  SECTION( "operator []" ) {
    seq[ 1 ] = 'G';
    REQUIRE( seq.sequence == "AGCT" );
  }

  SECTION( "operator == and !=" ) {
    REQUIRE( Sequence< DNA >( "ATCG" ) == Sequence< DNA >( "ATCG" ) );
    REQUIRE( Sequence< DNA >( "ATTG" ) != Sequence< DNA >( "ATCG" ) );

    REQUIRE( Sequence< DNA >( "NNTT" ) == Sequence< DNA >( "CGTT" ) );
    REQUIRE( Sequence< DNA >( "AUUT" ) == Sequence< DNA >( "ATTT" ) );

    REQUIRE( Sequence< DNA >( "ATT" ) != Sequence< DNA >( "ATTG" ) );
  }

  SECTION( "complement" ) {
    REQUIRE( Sequence< DNA >( "ACCT" ).Complement().sequence == "TGGA" );
    REQUIRE( Sequence< DNA >( "RGGN" ).Complement().sequence == "YCCN" );
  }

  SECTION( "reverse" ) {
    Sequence< DNA > rev = seq.Reverse();
    REQUIRE( seq.sequence == "ACCT" );
    REQUIRE( rev.sequence == "TCCA" );
  }
}
