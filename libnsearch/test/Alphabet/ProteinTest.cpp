#include <catch.hpp>

#include <nsearch/Alphabet/Protein.h>
const size_t BitMapPolicy< Protein >::NumBits;

TEST_CASE( "Protein" ) {
  SECTION( "Bitmapping" ) {
    REQUIRE( BitMapPolicy< Protein >::BitMap( 'C' ) >= 0 );

    // E is bitmapped the same as Q
    REQUIRE( BitMapPolicy< Protein >::BitMap( 'E' ) ==
             BitMapPolicy< Protein >::BitMap( 'Q' ) );

    REQUIRE( BitMapPolicy< Protein >::BitMap( 'X' ) == -1 );
  }

  SECTION( "Scoring and matching" ) {
    // BLOSUM62
    REQUIRE( ScorePolicy< Protein >::Score( 'H', 'H' ) == 8 );
    REQUIRE( ScorePolicy< Protein >::Score( 'P', 'T' ) == -1 );
    REQUIRE( ScorePolicy< Protein >::Score( 'B', 'N' ) == 3  );
    REQUIRE( ScorePolicy< Protein >::Score( 'C', 'K' ) == -3  );

    REQUIRE( MatchPolicy< Protein >::Match( 'A', 'A' ) == true );
    REQUIRE( MatchPolicy< Protein >::Match( 'Q', 'E' ) == false );
  }
}
