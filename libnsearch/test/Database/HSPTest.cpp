#include <catch.hpp>

#include <nsearch/Database/HSP.h>

#include <vector>

TEST_CASE( "HSP" )  {
  SECTION( "Basic methods" ) {
    HSP hsp( 5, 6, 11, 15, -62 );

    REQUIRE( hsp.Length() == 5 ); // 11-125
    REQUIRE( hsp.Score() == -62 );
  }

  SECTION( "operator <" ) {
    // Uses score
    HSP hsp1( 1, 1, 5, 5, 12 );
    HSP hsp2( 1, 1, 5, 5, 11 );

    REQUIRE( hsp2 < hsp1 );
  }

  SECTION( "IsOverlappinng" ) {
    HSP hsp( 1, 2, 25, 27 );

    REQUIRE( hsp.IsOverlapping( HSP( 1, 2, 55, 57 ) ) == true );
    REQUIRE( hsp.IsOverlapping( HSP( 2, 3, 55, 57 ) ) == true );
    REQUIRE( hsp.IsOverlapping( HSP( 3, 4, 55, 57 ) ) == false );
    REQUIRE( hsp.IsOverlapping( HSP( 3, 4, 20, 28 ) ) == true );
    REQUIRE( hsp.IsOverlapping( HSP( 3, 4, 20, 25 ) ) == true );
    REQUIRE( hsp.IsOverlapping( HSP( 3, 4, 20, 24 ) ) == false );
    REQUIRE( hsp.IsOverlapping( HSP( 0, 100, 0, 100 ) ) == true );
  }

  SECTION( "DistanceTo" ) {
    {
      HSP hsp1( 1, 1, 25, 25 );
      HSP hsp2( 2, 2, 26, 26 );
      REQUIRE( hsp1.DistanceTo( hsp2 ) == 0 );
    }
    {
      HSP hsp1( 1, 2, 25, 26 );
      HSP hsp2( 5, 10, 30, 35 );
      // x: 2 -> 5 = 2 empty cells
      // y: 26 -> 30 = 3 empty cells
      REQUIRE( hsp1.DistanceTo( hsp2 ) == ( size_t )sqrt( 2 * 2 + 3 * 3 )  );
      REQUIRE( hsp1.DistanceTo( hsp2 ) == hsp2.DistanceTo( hsp1 ) );
    }
  }
}
