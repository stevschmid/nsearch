#include <catch.hpp>

#include <nsearch/Alignment/BandedAlign.h>
#include <nsearch/Alphabet/DNA.h>
#include <nsearch/Sequence.h>

TEST_CASE( "BandedAlign" ) {
  Cigar cigar;

  SECTION( "Basic" ) {
    Sequence< DNA > a = "TATAATGTTTACATTGG";
    Sequence< DNA > b = "TATAATGACACTGG";

    BandedAlignParams  bap;
    BandedAlign< DNA > ba( bap );

    int score = ba.Align( a, b, &cigar );
    // TATAATGTTTACATTGG
    // |||||||   |||.|||
    // TATAATG---ACACTGG
    REQUIRE( cigar.ToString() == "7=3I3=1X3=" );
    const int matchScore = 2;
    const int mismatchScore = -4;
    REQUIRE( score == 13 * matchScore + 1 * bap.interiorGapOpenScore +
                        3 * bap.interiorGapExtendScore + 1 * mismatchScore );
  }

  SECTION( "Gap Penalties" ) {
    Sequence< DNA > a = "GGATCCTA";
    Sequence< DNA > b = "ATCGTA";

    {
      // Default: terminal gaps are not penalized heavily
      BandedAlign< DNA > ba;
      ba.Align( a, b, &cigar );
      REQUIRE( cigar.ToString() == "2I3=1X2=" );
    }

    {
      // Penalize terminal gaps heavily
      BandedAlignParams bap;
      bap.terminalGapOpenScore = bap.interiorGapOpenScore * 2;
      BandedAlign< DNA > ba( bap );
      ba.Align( a, b, &cigar );
      // GGATCCTA
      // A--TCGTA
      REQUIRE( cigar.ToString() == "1X2I2=1X2=" );
    }
  }

  SECTION( "Long tails" ) {
    BandedAlign< DNA > ba;
    ba.Align( "ATCGGGGGGGGGGGGGGGGGGGGGGG", "CGG", &cigar );
    REQUIRE( cigar.ToString() == "2I3=21I" );
    ba.Align( "GGGGTATAAAATTT", "TTTTTTTTGGGGTATAAAA", &cigar );
    REQUIRE( cigar.ToString() == "8D11=3I" );
  }

  SECTION( "Edge cases" ) {
    BandedAlign< DNA > ba;
    ba.Align( "", "", &cigar );
    REQUIRE( cigar.ToString() == "" );

    ba.Align( "A", "", &cigar );
    REQUIRE( cigar.ToString() == "1I" );

    ba.Align( "", "T", &cigar );
    REQUIRE( cigar.ToString() == "1D" );
  }

  SECTION( "Offsets" ) {
    BandedAlign< DNA > ba;
    SECTION( "Going forward" ) {
      ba.Align( "TTTTATCGGTAT", "GGCGGTAT", &cigar, AlignmentDirection::Forward, 0,
                0 );
      REQUIRE( cigar.ToString() == "4I2X6=" );

      ba.Align( "TTTTATCGGTAT", "GGCGGTAT", &cigar, AlignmentDirection::Forward, 4,
                2 );
      REQUIRE( cigar.ToString() == "2I6=" );
    }

    SECTION( "Going backwards" ) {
      // GGATGA-
      // --ATGAA
      ba.Align( "GGATGA", "ATGAA", &cigar, AlignmentDirection::Reverse, 6, 5 );
      REQUIRE( cigar.ToString() == "2I4=1D" );

      // GGATGA
      // --ATG-
      ba.Align( "GGATGA", "ATGAA", &cigar, AlignmentDirection::Reverse, 6, 3 );
      REQUIRE( cigar.ToString() == "2I3=1I" );
    }
  }

  // Breaking cases
  SECTION( "Breaking case when first row is not initialized properly (beyond "
           "bandwidth)" ) {
    Sequence< DNA >    A = "AAAAAAAAAAAAAAA";
    Sequence< DNA >    B = "CCCCCCAAAAAAAAA";
    BandedAlign< DNA > ba;

    ba.Align( A, B, &cigar );
    REQUIRE( cigar.ToString() == "6D9=6I" );

    A = "CCCCCCCCCCCCCCC";
    B = "CCCCCCAAAAAAAAA";
    ba.Align( A, B, &cigar );
    REQUIRE( cigar.ToString() == "9I6=9D" );
  }

  SECTION( "Breaking case when startA >> lenA" ) {
    Sequence< DNA >    A = "ATGCC";
    Sequence< DNA >    B = "TTTATGCC";
    BandedAlign< DNA > ba;

    ba.Align( A, B, &cigar, AlignmentDirection::Forward, 6, 3 );
    REQUIRE( cigar.ToString() == "5D" );
  }

  SECTION( "Breaking case: Improper reset of vertical gap at 0,0" ) {
    Cigar              cigar1, cigar2;
    BandedAlign< DNA > ba;

    // Align first with fresh alignment cache
    int score1 =
      ba.Align( "ATGCC", "TTTTAGCC", &cigar1, AlignmentDirection::Forward, 1, 1 );
    REQUIRE( cigar1.ToString() == "1=3X3D" );

    // This alignment will set mVerticalGaps[0] to a low value, which will be
    // extended upon subsequently if we don't reset
    ba.Align( "ATGCC", "A", &cigar, AlignmentDirection::Forward );

    // Test with the "leaky" vgap
    int score2 =
      ba.Align( "ATGCC", "TTTTAGCC", &cigar2, AlignmentDirection::Forward, 1, 1 );

    REQUIRE( cigar1.ToString() == cigar2.ToString() );
    REQUIRE( score1 == score2 );
  }
}
