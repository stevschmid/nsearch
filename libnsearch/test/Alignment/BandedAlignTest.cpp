#include <catch.hpp>

#include <nsearch/Alignment/BandedAlign.h>

TEST_CASE( "BandedAlign" )  {
  Cigar cigar;

  SECTION( "Basic" ) {
    Sequence a = "TATAATGTTTACATTGG";
    Sequence b = "TATAATGACACTGG";

    BandedAlignParams bap;
    BandedAlign ba( bap );
    int score = ba.Align( a, b, &cigar );
    // TATAATGTTTACATTGG
    // |||||||   |||.|||
    // TATAATG---ACACTGG
    REQUIRE( cigar.ToString() == "7M3I3M1X3M" );
    REQUIRE( score ==
        13 * bap.matchScore +
        1 * bap.interiorGapOpenScore +
        3 * bap.interiorGapExtendScore +
        1 * bap.mismatchScore );
  }

  SECTION( "Gap Penalties" ) {
    Sequence a = "GGATCCTA";
    Sequence b = "ATCGTA";

    {
      // Default: terminal gaps are not penalized heavily
      BandedAlign ba;
      ba.Align( a, b, &cigar );
      REQUIRE( cigar.ToString() == "2I3M1X2M" );
    }

    {
      // Penalize terminal gaps heavily
      BandedAlignParams bap;
      bap.terminalGapOpenScore = bap.interiorGapOpenScore * 2;
      BandedAlign ba( bap );
      ba.Align( a, b, &cigar );
      // GGATCCTA
      // A--TCGTA
      REQUIRE( cigar.ToString() == "1X2I2M1X2M" );
    }
  }

  SECTION( "Long tails" ) {
    BandedAlign ba;
    ba.Align( "ATCGGGGGGGGGGGGGGGGGGGGGGG", "CGG", &cigar );
    REQUIRE( cigar.ToString() == "2I3M21I" );
    ba.Align( "GGGGTATAAAATTT", "TTTTTTTTGGGGTATAAAA", &cigar );
    REQUIRE( cigar.ToString() == "8D11M3I" );
  }

  SECTION( "Edge cases" ) {
    BandedAlign ba;
    ba.Align( "", "", &cigar );
    REQUIRE( cigar.ToString() == "" );

    ba.Align( "A", "", &cigar );
    REQUIRE( cigar.ToString() == "1I" );

    ba.Align( "", "T", &cigar );
    REQUIRE( cigar.ToString() == "1D" );
  }

  SECTION( "Offsets" ) {
    BandedAlign ba;
    SECTION( "Going forward") {
      ba.Align( "TTTTATCGGTAT", "GGCGGTAT", &cigar, AlignmentDirection::fwd, 0, 0 );
      REQUIRE( cigar.ToString() == "4I2X6M" );

      ba.Align( "TTTTATCGGTAT", "GGCGGTAT", &cigar, AlignmentDirection::fwd, 4, 2 );
      REQUIRE( cigar.ToString() == "2I6M" );
    }

    SECTION( "Going backwards" ) {
      // GGATGA-
      // --ATGAA
      ba.Align( "GGATGA", "ATGAA", &cigar, AlignmentDirection::rev, 6, 5 );
      REQUIRE( cigar.ToString() == "2I4M1D" );

      // GGATGA
      // --ATG-
      ba.Align( "GGATGA", "ATGAA", &cigar, AlignmentDirection::rev, 6, 3 );
      REQUIRE( cigar.ToString() == "2I3M1I" );
    }
  }

  // Breaking cases
  SECTION( "Breaking case when first row is not initialized properly (beyond bandwidth)" ) {
    Sequence A = "AAAAAAAAAAAAAAA";
    Sequence B = "CCCCCCAAAAAAAAA";
    BandedAlign ba;
    ba.Align( A, B, &cigar );
    REQUIRE( cigar.ToString() == "6D9M6I" );

    A = "CCCCCCCCCCCCCCC";
    B = "CCCCCCAAAAAAAAA";
    ba.Align( A, B, &cigar );
    REQUIRE( cigar.ToString() == "9I6M9D" );
  }

  SECTION( "Breaking case when startA >> lenA" ) {
    Sequence A = "ATGCC";
    Sequence B = "TTTATGCC";
    BandedAlign ba;
    ba.Align( A, B, &cigar, AlignmentDirection::fwd, 6, 3 );
    REQUIRE( cigar.ToString() == "5D" );
  }

  SECTION( "Breaking case: Improper reset of vertical gap at 0,0" ) {
    Cigar cigar1, cigar2;
    BandedAlign ba;

    // Align first with fresh alignment cache
    int score1 = ba.Align( "ATGCC", "TTTTAGCC", &cigar1, AlignmentDirection::fwd, 1, 1 );
    REQUIRE( cigar1.ToString() == "1M3X3D" );

    // This alignment will set mVerticalGaps[0] to a low value, which will be extended
    // upon subsequently if we don't reset
    ba.Align( "ATGCC", "A", &cigar,  AlignmentDirection::fwd );

    // Test with the "leaky" vgap
    int score2 = ba.Align( "ATGCC", "TTTTAGCC", &cigar2,  AlignmentDirection::fwd, 1, 1 );

    REQUIRE( cigar1.ToString() == cigar2.ToString() );
    REQUIRE( score1 == score2 );
  }
}
