#include <catch.hpp>

#include <nsearch/Aligner.h>

TEST_CASE( "Aligner" )  {
  Aligner aligner( 2, -3, 5, 2 );

  SECTION( "Global" ) {
    Sequence query( "GACTTAC");
    Sequence target( "CGTGAATTCAT" );

    SECTION( "Score" ) {
      REQUIRE( aligner.GlobalAlign( query, target ) == -14 );
    }

    SECTION( "Symmetry" ) {
      REQUIRE( aligner.GlobalAlign( query, target ) == aligner.GlobalAlign( target, query ) );
    }

    SECTION( "Alignment" ) {
      Alignment aln;
      aligner.GlobalAlign( query, target, &aln );
      REQUIRE( aln.cigarString() == "3D5M1D2M" );
      REQUIRE( aln.score == -14 );
      REQUIRE( aln.targetPos == 0 );
      REQUIRE( aln.queryPos == 0 );
    }
  }

  SECTION( "Local" ) {
    Sequence query( "TACGGGCCCGCTAC" );
    Sequence target( "TAGCCCTATCGGTCA" );
    Alignment aln;

    SECTION( "Score" ) {
      Aligner aligner( 5, -4, 4, 1 );
      REQUIRE( aligner.LocalAlign( query, target, &aln ) == 27 );
      REQUIRE( aln.cigarString() == "2M3I4M2I2M7S" );
    }

    SECTION( "Symmetry" ) {
      REQUIRE( aligner.LocalAlign( query, target ) == aligner.LocalAlign( target, query ) );
    }

    SECTION( "Gaps" ) {
      REQUIRE( Aligner( 5, -4, 4, 1 ).LocalAlign( query, target ) == 27 );
      REQUIRE( Aligner( 5, -4, 10, 1 ).LocalAlign( query, target ) == 20 );
    }

    SECTION( "Ambiguous Nucleotides" ) {
      Sequence query = "ATGCTGGTACCTGGGAT";
      Sequence target =       "TGGYAWNNV"; // W matches A or T, here C (mismatch)
      REQUIRE( aligner.LocalAlign( query, target ) == 8 * 2 - 1 * 3 );
    }

    SECTION( "Alignment" ) {
      Sequence query =         "TACCTGG";
      Sequence target = "ATGCTGGTACCTGGGAT";
      REQUIRE( aligner.LocalAlign( query, target, &aln ) == 14 );
      REQUIRE( aln.score == 14 );
      REQUIRE( aln.targetPos == 7 );
      REQUIRE( aln.queryPos == 0 );

      REQUIRE( aligner.LocalAlign( target, query, &aln ) == 14 );
      REQUIRE( aln.targetPos == 0 );
      REQUIRE( aln.queryPos == 7 );
    }

    SECTION( "Clipping" ) {
      Sequence query =      "TGGT";
      Sequence target = "ATGCTGGTACCTGG";
      REQUIRE( aligner.LocalAlign( query, target, &aln ) == 8 );
    }
  }
}
