#include <catch.hpp>

#include <nsearch/Alignment/Ranges.h>

TEST_CASE( "Ranges" )  {
  Ranges r;

  SECTION( "Simple overlap" ) {
    r.Add( 0, 10 );
    r.Add( 5, 15 );
    REQUIRE( r.NumRanges() == 1 );
    REQUIRE( r.Size() == 16 );
  }

  SECTION( "Overlapping" ) {
    r.Add( 1, 3 );
    r.Add( 5, 7 );
    r.Add( 4, 10 );
    REQUIRE( r.NumRanges() == 2 );
    REQUIRE( r.Size() == 10 );

    r.Add( 15, 19 );
    REQUIRE( r.NumRanges() == 3 );
    REQUIRE( r.Size() == 15 );

    r.Add( 1, 20 );
    REQUIRE( r.NumRanges() == 1 );
    REQUIRE( r.Size() == 20 );
  }

  SECTION( "operator<" ) {
    Ranges small;
    small.Add( 1, 3 ); // 3
    Ranges medium;
    medium.Add( 2, 5 ); // 4
    Ranges big;
    big.Add( 3, 9 ); // 7

    REQUIRE( small < medium );
    REQUIRE( medium < big );
    REQUIRE( small < big );
  }
}
