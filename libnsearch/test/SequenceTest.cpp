#include <catch.hpp>

#include <nsearch/Sequence.h>

TEST_CASE( "Bla" )  {
  Sequence seq( "@Id", "ACCT", "JJ::" );

  REQUIRE( seq.Length() == 4 );
  REQUIRE( seq.identifier == "@Id" );
  REQUIRE( seq.sequence == "ACCT" );
  REQUIRE( seq.quality == "JJ::" );

  SECTION( "subsequence" ) {
    Sequence sub = seq.Subsequence( 2, 2 );
    REQUIRE( sub.sequence == "CT" );
    REQUIRE( sub.quality == "::" );

    Sequence sub2 = seq.Subsequence( 1 );
    REQUIRE( sub2.sequence == "CCT" );
  }

  SECTION( "operator +" ) {
    Sequence m = seq + Sequence( "@Id2", "GCTT",  "BBAJ" );

    REQUIRE( m.identifier == "@Id" );
    REQUIRE( m.sequence == "ACCTGCTT" );
    REQUIRE( m.quality == "JJ::BBAJ" );
  }

  SECTION( "operator []" ) {
    seq[ 1 ] = 'G';
    REQUIRE( seq.sequence == "AGCT" );
  }

  SECTION( "complement" ) {
    Sequence com = seq.Complement();
    REQUIRE( seq.sequence == "ACCT" );
    REQUIRE( com.sequence == "TGGA" );
  }

  SECTION( "reverse" ) {
    Sequence rev = seq.Reverse();
    REQUIRE( seq.sequence == "ACCT" );
    REQUIRE( rev.sequence == "TCCA" );
  }

  SECTION( "reverse complement" ) {
    Sequence rco = seq.ReverseComplement();
    REQUIRE( seq.sequence == "ACCT" );
    REQUIRE( rco.sequence == "AGGT" );
  }
}
