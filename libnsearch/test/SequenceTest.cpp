#include <catch.hpp>

#include <nsearch/Sequence.h>

TEST_CASE( "Sequence" ) {
  Sequence seq( "Id", "ACCT", "JJ::" );

  REQUIRE( seq.Length() == 4 );
  REQUIRE( seq.identifier == "Id" );
  REQUIRE( seq.sequence == "ACCT" );
  REQUIRE( seq.quality == "JJ::" );

  SECTION( "subsequence" ) {
    Sequence sub = seq.Subsequence( 2, 2 );
    REQUIRE( sub.sequence == "CT" );
    REQUIRE( sub.quality == "::" );

    REQUIRE( seq.Subsequence( 2 ).sequence == "CT" );
    REQUIRE( seq.Subsequence( 4 ).sequence == "" );

    Sequence sub3 = Sequence( "ATTG" ).Subsequence( 10, 10 );
    REQUIRE( sub3.sequence == "" );
  }

  SECTION( "operator +" ) {
    Sequence m = seq + Sequence( "Id2", "GCTT", "BBAJ" );

    REQUIRE( m.identifier == "Id" );
    REQUIRE( m.sequence == "ACCTGCTT" );
    REQUIRE( m.quality == "JJ::BBAJ" );
  }

  SECTION( "operator []" ) {
    seq[ 1 ] = 'G';
    REQUIRE( seq.sequence == "AGCT" );
  }

  SECTION( "operator == and !=" ) {
    REQUIRE( Sequence( "ATCG" ) == Sequence( "ATCG" ) );
    REQUIRE( Sequence( "ATTG" ) != Sequence( "ATCG" ) );

    REQUIRE( Sequence( "NNTT" ) == Sequence( "CGTT" ) );
    REQUIRE( Sequence( "AUUT" ) == Sequence( "ATTT" ) );

    REQUIRE( Sequence( "ATT" ) != Sequence( "ATTG" ) );
  }

  SECTION( "complement" ) {
    REQUIRE( Sequence( "ACCT" ).Complement().sequence == "TGGA" );
    REQUIRE( Sequence( "RGGN" ).Complement().sequence == "YCCN" );
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
