#include <catch.hpp>

#include <nsearch/Alignment/ExtendAlign.h>
#include <nsearch/Alphabet/DNA.h>
#include <nsearch/Sequence.h>

TEST_CASE( "ExtendAlign" ) {
  size_t bestA, bestB;
  Cigar  cigar;
  int    score;

  ExtendAlign< DNA > ea;

  SECTION( "Gaps" ) {
    Sequence< DNA > a = "GATTGCGGGG";
    Sequence< DNA > b = "GAGCGGT";

    score =
      ea.Extend( a, b, &bestA, &bestB, &cigar, AlignmentDirection::Forward, 0, 0 );
    REQUIRE( cigar.ToString() == "2M" );

    ExtendAlignParams eap;
    eap.gapOpenScore = eap.gapExtendScore - 1;

    ea = ExtendAlign< DNA >( eap );
    score =
      ea.Extend( a, b, &bestA, &bestB, &cigar, AlignmentDirection::Forward, 0, 0 );
    REQUIRE( cigar.ToString() == "2M2I4M" );
  }

  SECTION( "Forwards extend" ) {
    Sequence< DNA > a = "ATCGG";
    Sequence< DNA > b = "ATCGT";

    score =
      ea.Extend( a, b, &bestA, &bestB, &cigar, AlignmentDirection::Forward, 0, 0 );
    REQUIRE( bestA == 3 );
    REQUIRE( bestB == 3 );
    REQUIRE( cigar.ToString() == "4M" );

    score =
      ea.Extend( a, b, &bestA, &bestB, &cigar, AlignmentDirection::Forward, 3, 3 );
    REQUIRE( bestA == 3 );
    REQUIRE( bestB == 3 );
    REQUIRE( cigar.ToString() == "1M" );

    score =
      ea.Extend( a, b, &bestA, &bestB, &cigar, AlignmentDirection::Forward, 4, 4 );
    REQUIRE( bestA == 4 );
    REQUIRE( bestB == 4 );
    REQUIRE( cigar.ToString() == "" );
  }

  SECTION( "Backwards extend" ) {
    Sequence< DNA > a = "ATCGGTTG";
    Sequence< DNA > b = "TCGGTAT";

    score =
      ea.Extend( a, b, &bestA, &bestB, &cigar, AlignmentDirection::Reverse, 3, 2 );
    REQUIRE( bestA == 1 );
    REQUIRE( bestB == 0 );

    score =
      ea.Extend( a, b, &bestA, &bestB, &cigar, AlignmentDirection::Reverse, 1, 0 );
    REQUIRE( bestA == 1 );
    REQUIRE( bestB == 0 );
  }
}
