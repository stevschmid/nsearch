#include <catch.hpp>

#include <nsearch/Alignment/Cigar.h>

TEST_CASE( "CigarTest" ) {
  SECTION( "Constructing" ) {
    Cigar cigar1;
    Cigar cigar2( "5M10I" );
    REQUIRE( cigar2[ 0 ] == CigarEntry( 5, CigarOp::Match ) );
    REQUIRE( cigar2[ 1 ] == CigarEntry( 10, CigarOp::Insertion ) );
  }

  SECTION( "Adding" ) {
    Cigar cigar;
    cigar.Add( { 2, CigarOp::Match } );
    cigar.Add( { 5, CigarOp::Match } );
    cigar.Add( { 20, CigarOp::Deletion } );
    cigar.Add( { 3, CigarOp::Mismatch } );
    cigar.Add( { 1, CigarOp::Match } );
    cigar.Add( CigarOp::Match );
    cigar.Add( CigarOp::Match );

    REQUIRE( cigar.ToString() == "7M20D3X3M" );
  }

  SECTION( "Identity" ) {
    REQUIRE( Cigar( "50I14M2X4M25D" ).Identity() == float( 18 ) / float( 20 ) );
    REQUIRE( Cigar( "2M" ).Identity() == 1.0f );
    REQUIRE( Cigar( "1X1M" ).Identity() == 0.5f );
  }
}
