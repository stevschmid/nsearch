#include <catch.hpp>

#include <nsearch/Alignment/HitTracker.h>

TEST_CASE( "HitTracker" )  {
  HitTracker ht;

  SECTION( "Same diagonal joined" ) {

    ht.AddHit( 0, 10, 5 ); // diag1
    ht.AddHit( 0, 11, 3 ); // diag2
    ht.AddHit( 1, 11, 6 ); // diag1

    SegmentPairList sps = ht.List();
    REQUIRE( sps.size() == 2 );
    REQUIRE( std::find( sps.begin(), sps.end(), SegmentPair( 0, 10, 7 ) ) != sps.end() );
    REQUIRE( std::find( sps.begin(), sps.end(), SegmentPair( 0, 11, 3 ) ) != sps.end() );
  }

  SECTION( "Same diagonal yet disjoint" ) {
    ht.AddHit( 10, 5, 5 ); // diag1
    ht.AddHit( 14, 9, 2 ); // diag1 joined

    ht.AddHit( 20, 15, 2 ); // diag1 disjoint

    ht.AddHit( 3, 18, 1 );
    ht.AddHit( 55, 2, 10 );

    SegmentPairList sps = ht.List();
    std::sort( sps.begin(), sps.end(), []( const SegmentPair &l, const SegmentPair &r ) {
      return l.s1 < r.s1;
    });

    REQUIRE( sps.size() == 4 );
    REQUIRE( sps[ 0 ] == SegmentPair( 3, 18, 1 ) );
    REQUIRE( sps[ 1 ] == SegmentPair( 10, 5, 6 ) );
    REQUIRE( sps[ 2 ] == SegmentPair( 20, 15, 2 ) );
    REQUIRE( sps[ 3 ] == SegmentPair( 55, 2, 10 ) );
  }
}
