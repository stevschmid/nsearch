#include <catch.hpp>

#include <nsearch/Alignment/RangeMerger.h>

TEST_CASE( "RangeMerger" )  {
  RangeMerger rm;

  SECTION( "Simple overlap" ) {
    rm.Add( 0, 10 );
    rm.Add( 5, 15 );
    REQUIRE( rm.NumRanges() == 1 );
    REQUIRE( rm.Size() == 16 );
  }

  SECTION( "Overlapping" ) {
    rm.Add( 1, 3 );
    rm.Add( 5, 7 );
    rm.Add( 4, 10 );
    REQUIRE( rm.NumRanges() == 2 );
    REQUIRE( rm.Size() == 10 );

    rm.Add( 15, 19 );
    REQUIRE( rm.NumRanges() == 3 );
    REQUIRE( rm.Size() == 15 );

    rm.Add( 1, 20 );
    REQUIRE( rm.NumRanges() == 1 );
    REQUIRE( rm.Size() == 20 );
  }

  SECTION( "operator<" ) {
    RangeMerger small;
    small.Add( 1, 3 ); // 3
    RangeMerger medium;
    medium.Add( 2, 5 ); // 4
    RangeMerger big;
    big.Add( 3, 9 ); // 7

    REQUIRE( small < medium );
    REQUIRE( medium < big );
    REQUIRE( small < big );
  }
}
