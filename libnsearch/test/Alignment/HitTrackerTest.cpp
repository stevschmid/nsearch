#include <catch.hpp>

#include <nsearch/Alignment/HitTracker.h>

TEST_CASE( "HitTracker" )  {
  HitTracker ht;

  SECTION( "Same diagonal joined" ) {

    ht.AddHit( 0, 10, 5 ); // diag1
    ht.AddHit( 0, 11, 3 ); // diag2
    ht.AddHit( 1, 11, 6 ); // diag1

    REQUIRE( ht.Seeds().size() == 2 );
    REQUIRE( ht.Seeds()[ 0 ] == Seed( 0, 10, 7 ) );
    REQUIRE( ht.Seeds()[ 1 ] == Seed( 0, 11, 3 ) );
  }

  SECTION( "Same diagonal yet disjoint" ) {
    ht.AddHit( 10, 5, 5 ); // diag1
    ht.AddHit( 14, 9, 2 ); // diag1 joined

    ht.AddHit( 20, 15, 2 ); // diag1 disjoint

    ht.AddHit( 3, 18, 1 );
    ht.AddHit( 55, 2, 10 );

    SeedList seeds = ht.Seeds();
    std::sort( seeds.begin(), seeds.end(), []( const Seed &l, const Seed &r ) {
      return l.s1 < r.s1;
    });

    REQUIRE( seeds.size() == 4 );
    REQUIRE( seeds[ 0 ] == Seed( 3, 18, 1 ) );
    REQUIRE( seeds[ 1 ] == Seed( 10, 5, 6 ) );
    REQUIRE( seeds[ 2 ] == Seed( 20, 15, 2 ) );
    REQUIRE( seeds[ 3 ] == Seed( 55, 2, 10 ) );
  }
}
