#include <catch.hpp>

#include <nsearch/Alphabet/DNA.h>
const size_t BitMapPolicy< DNA >::NumBits;

TEST_CASE( "DNA" ) {
  SECTION( "Bitmapping" ) {
    REQUIRE( BitMapPolicy< DNA >::BitMap( 'U' ) ==
             BitMapPolicy< DNA >::BitMap( 'T' ) );
    REQUIRE( BitMapPolicy< DNA >::BitMap( 'N' ) == -1 );
    REQUIRE( BitMapPolicy< DNA >::NumBits == 2 );
  }

  SECTION( "Complement" ) {
    REQUIRE( ComplementPolicy< DNA >::Complement( 'A' ) == 'T' );
    REQUIRE( ComplementPolicy< DNA >::Complement( 'M' ) == 'K' );
  }

  SECTION( "Scoring and matching" ) {
    REQUIRE( ScorePolicy< DNA >::Score( 'A', 'A' ) == 2 );
    REQUIRE( ScorePolicy< DNA >::Score( 'A', 'N' ) == 2 );
    REQUIRE( ScorePolicy< DNA >::Score( 'G', 'C' ) == -4 );

    REQUIRE( MatchPolicy< DNA >::Match( 'A', 'A' ) == true );
    REQUIRE( MatchPolicy< DNA >::Match( 'A', 'T' ) == false );
    REQUIRE( MatchPolicy< DNA >::Match( 'A', 'N' ) == true );
    REQUIRE( MatchPolicy< DNA >::Match( 'R', 'C' ) == false );
    REQUIRE( MatchPolicy< DNA >::Match( 'R', 'G' ) == true );
  }
}
