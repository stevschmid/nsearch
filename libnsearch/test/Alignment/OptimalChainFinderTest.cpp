#include <catch.hpp>

#include <nsearch/Alignment/OptimalChainFinder.h>

TEST_CASE( "OptimalChainFinder" )  {
  {
    SeedList seeds;
    seeds.push_back( Seed( 3, 908, 15 ) );
    seeds.push_back( Seed( 19, 924, 21 ) );
    OptimalChainFinder ch( seeds );

    SeedList optimalChain = ch.OptimalChain();
    REQUIRE( optimalChain.size() == 2 );
    REQUIRE( optimalChain[ 0 ] == Seed( 3, 908, 15 ) );
    REQUIRE( optimalChain[ 1 ] == Seed( 19, 924, 21 ) );
  }

  {
    SeedList seeds;
    seeds.push_back( Seed( 0, 3, 1 ) );
    seeds.push_back( Seed( 1, 4, 5 ) );
    OptimalChainFinder ch( seeds );

    SeedList optimalChain = ch.OptimalChain();
    REQUIRE( optimalChain.size() == 2 );
    REQUIRE( optimalChain[ 0 ] == Seed( 0, 3, 1 ) );
    REQUIRE( optimalChain[ 1 ] == Seed( 1, 4, 5 ) );
  }

  {
    SeedList seeds;
    seeds.push_back( Seed( 0, 0, 3 ) );
    seeds.push_back( Seed( 3, 0, 4 ) );
    seeds.push_back( Seed( 7, 4, 3 ) );
    OptimalChainFinder ch( seeds );

    SeedList optimalChain = ch.OptimalChain();
    std::cout << optimalChain.size() << std::endl;
    REQUIRE( optimalChain.size() == 2 );
    REQUIRE( optimalChain[ 0 ] == Seed( 3, 0, 4 ) );
    REQUIRE( optimalChain[ 1 ] == Seed( 7, 4, 3 ) );
  }

  SECTION( "Seed order should not matter" ) {
    SeedList seeds;
    seeds.push_back( Seed( 0, 0, 2 ) );
    seeds.push_back( Seed( 2, 2, 2 ) );
    OptimalChainFinder ch( seeds );
    REQUIRE( ch.OptimalChain().size() == 2 );

    std::reverse( seeds.begin(), seeds.end() );
    OptimalChainFinder chRev( seeds );
    REQUIRE( chRev.OptimalChain().size() == 2 );
  }

  {
    // https://github.com/seqan/seqan/issues/2082
    SeedList seeds;
    seeds.push_back( Seed( 0, 0, 14 ) );
    seeds.push_back( Seed( 13, 20, 13 ) );
    OptimalChainFinder ch( seeds );

    SeedList optimalChain = ch.OptimalChain();
    REQUIRE( optimalChain.size() == 1 );
    REQUIRE( optimalChain[ 0 ] == Seed( 0, 0, 14 ) );
  }

  {
    SeedList seeds;
    seeds.push_back( Seed( 0, 0, 1 ) );
    seeds.push_back( Seed( 2, 4, 2 ) );
    seeds.push_back( Seed( 2, 0, 3 ) );
    seeds.push_back( Seed( 7, 2, 6 ) );

    seeds.push_back( Seed( 3, 14, 3 ) );
    seeds.push_back( Seed( 6, 14, 4 ) );
    seeds.push_back( Seed( 10, 18, 3 ) );

    seeds.push_back( Seed( 15, 11, 4 ) );
    seeds.push_back( Seed( 16, 5, 2 ) );
    seeds.push_back( Seed( 14, 21, 2 ) );

    seeds.push_back( Seed( 8, 22, 6 ) );
    OptimalChainFinder ch( seeds );

    SeedList optimalChain = ch.OptimalChain();
    REQUIRE( optimalChain.size() == 5 );
    REQUIRE( optimalChain[ 0 ] == Seed( 0, 0, 1 ) );
    REQUIRE( optimalChain[ 1 ] == Seed( 2, 4, 2 ) );
    REQUIRE( optimalChain[ 2 ] == Seed( 6, 14, 4 ) );
    REQUIRE( optimalChain[ 3 ] == Seed( 10, 18, 3 ) );
    REQUIRE( optimalChain[ 4 ] == Seed( 14, 21, 2 ) );
  }
}
