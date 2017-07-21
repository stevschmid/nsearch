#include <catch.hpp>

#include <nsearch/Alignment/ExtendAlign.h>

TEST_CASE( "ExtendAlign" )  {
  size_t bestA, bestB;
  Cigar cigar;
  int score;

  ExtendAlign ea;

  SECTION( "Gaps" ) {
    Sequence a = "GATTGCGGGG";
    Sequence b = "GAGCGGT";

    score = ea.Extend( a, b, &bestA, &bestB, &cigar,
        AlignmentDirection::forwards,
        0, 0 );
    REQUIRE( cigar.ToString() == "2M" );

    ExtendAlignParams eap;
    eap.gapOpenScore = eap.gapExtendScore -1;
    ea = ExtendAlign( eap );
    score = ea.Extend( a, b, &bestA, &bestB, &cigar,
        AlignmentDirection::forwards,
        0, 0 );
    REQUIRE( cigar.ToString() == "2M2I4M" );
  }

  SECTION( "Forwards extend" ) {
    Sequence a = "ATCGG";
    Sequence b = "ATCGT";

    score = ea.Extend( a, b, &bestA, &bestB, &cigar,
        AlignmentDirection::forwards,
        0, 0 );
    REQUIRE( bestA == 3 );
    REQUIRE( bestB == 3 );
    REQUIRE( cigar.ToString() == "4M" );

    score = ea.Extend( a, b, &bestA, &bestB, &cigar,
        AlignmentDirection::forwards,
        3, 3 );
    REQUIRE( bestA == 3 );
    REQUIRE( bestB == 3 );
    REQUIRE( cigar.ToString() == "1M" );

    score = ea.Extend( a, b, &bestA, &bestB, &cigar,
        AlignmentDirection::forwards,
        4, 4 );
    REQUIRE( bestA == 4 );
    REQUIRE( bestB == 4 );
    REQUIRE( cigar.ToString() == "" );
  }

  SECTION( "Backwards extend" ) {
    Sequence a = "ATCGGTTG";
    Sequence b = "TCGGTAT";

    score = ea.Extend( a, b, &bestA, &bestB, &cigar,
        AlignmentDirection::backwards,
        3, 2 );
    REQUIRE( bestA == 1 );
    REQUIRE( bestB == 0 );

    score = ea.Extend( a, b, &bestA, &bestB, &cigar,
        AlignmentDirection::backwards,
        1, 0 );
    REQUIRE( bestA == 1 );
    REQUIRE( bestB == 0 );
  }

}
