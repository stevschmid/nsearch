#include <catch.hpp>

#include <nsearch/Alignment/Cigar.h>

TEST_CASE( "CigarTest" ) {
  SECTION( "Constructing" ) {
    Cigar cigar1;
    Cigar cigar2( "5=10I" );
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

    REQUIRE( cigar.ToString() == "7=20D3X3=" );
  }

  SECTION( "Identity" ) {
    REQUIRE( Cigar( "50I14=2X4=25D" ).Identity() == float( 18 ) / float( 20 ) );
    REQUIRE( Cigar( "2=" ).Identity() == 1.0f );
    REQUIRE( Cigar( "1X1=" ).Identity() == 0.5f );
  }
}
