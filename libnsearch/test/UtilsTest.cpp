#include <catch.hpp>

#include <nsearch/Utils.h>

TEST_CASE( "Utils" )  {
  SECTION( "IsBlank" ) {
    REQUIRE( IsBlank( "" ) == true );
    REQUIRE( IsBlank( " " ) == true );
    REQUIRE( IsBlank( "\r\n" ) == true );
    REQUIRE( IsBlank( "!" ) == false );
    REQUIRE( IsBlank( "\r\nA" ) == false );
  }

  SECTION( "UpcaseString" ) {
    std::string str = "AcGt";
    UpcaseString( str );
    REQUIRE( str == "ACGT" );
  }

  SECTION( "DoNucleotidesMatch" ) {
    REQUIRE( DoNucleotidesMatch( 'A', 'A' ) == true );
    REQUIRE( DoNucleotidesMatch( 'C', 'T' ) == false );
    REQUIRE( DoNucleotidesMatch( 'G', 'R' ) == true );
    REQUIRE( DoNucleotidesMatch( 'A', 'Y' ) == false );
    REQUIRE( DoNucleotidesMatch( 'M', 'Y' ) == true );
  }

  SECTION( "NucleotideComplement" ) {
    REQUIRE( NucleotideComplement( 'A' ) == 'T' );
    REQUIRE( NucleotideComplement( 'U' ) == 'A' );
    REQUIRE( NucleotideComplement( 'R' ) == 'Y' );
    REQUIRE( NucleotideComplement( 'N' ) == 'N' );
  }

  SECTION( "Coverage" ) {
    SECTION( "Fraction" ) {
      Coverage cov( 50 );

      REQUIRE( cov.CoveredFraction() == 0.0f );

      cov.Add( 1, 10 );
      cov.Add( 5, 15 );
      REQUIRE( cov.NumNonOverlaps() == 1 );
      REQUIRE( cov.CoveredSize() == 15 );
      REQUIRE( cov.CoveredFraction() == 0.30f );
    }

    SECTION( "Overlapping" ) {
      Coverage cov( 20 );

      cov.Add( 1, 3 );
      cov.Add( 5, 7 );
      cov.Add( 4, 10 );
      REQUIRE( cov.NumNonOverlaps() == 2 );
      REQUIRE( cov.CoveredSize() == 10 );

      cov.Add( 15, 19 );
      REQUIRE( cov.NumNonOverlaps() == 3 );
      REQUIRE( cov.CoveredSize() == 15 );

      cov.Add( 1, 20 );
      REQUIRE( cov.NumNonOverlaps() == 1 );
      REQUIRE( cov.CoveredSize() == 20 );
      REQUIRE( cov.CoveredFraction() == 1.0 );
    }

    SECTION( "operator<" ) {
      Coverage small( 10 );
      small.Add( 1, 3 ); // 3
      Coverage medium( 10 );
      medium.Add( 2, 5 ); // 4
      Coverage big( 10 );
      big.Add( 3, 9 ); // 7

      REQUIRE( small < medium );
      REQUIRE( medium < big );
      REQUIRE( small < big );

      std::set< Coverage > covs;
      covs.insert( medium );
      covs.insert( big );
      covs.insert( small );

      REQUIRE( (*covs.begin()).CoveredSize() == 3 );
      REQUIRE( (*covs.rbegin()).CoveredSize() == 7 );
    }
  }
}
