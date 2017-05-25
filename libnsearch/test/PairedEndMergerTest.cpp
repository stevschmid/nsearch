#include <catch.hpp>

#include <nsearch/PairedEndMerger.h>

TEST_CASE( "PairedEndMerger" )  {
  bool res;
  Sequence merged;

  SECTION( "non-staggered overlap" ) {
    PairedEndMerger merger( 5, 1.0 );
    Sequence fwd1 = Sequence( "fwd1", "ACTGGATGGA", "JJJJJJJJJJ" );
    Sequence rev1 = Sequence( "rev1",      "ATGGAATCCC", "JJJJJJJJJJ" ).ReverseComplement();

    res = merger.Merge( merged, fwd1, rev1 );
    REQUIRE( res == true );
    REQUIRE( merged == Sequence( "ACTGGATGGAATCCC" ) );
  }

  SECTION( "staggered overlap (merged sequence is trimmed)" ) {
    PairedEndMerger merger( 5, 1.0 );
    Sequence fwd1 = Sequence( "fwd1",      "ATCCCGGA", "JJJJJJJJJJ" );
    Sequence rev1 = Sequence( "rev1", "ATGGAATCCC", "JJJJJJJJJJ" ).ReverseComplement();

    res = merger.Merge( merged, fwd1, rev1 );
    REQUIRE( res == true );
    REQUIRE( merged == Sequence( "ATCCC" ) );
  }

  SECTION( "min bases requirement" ) {
    PairedEndMerger merger( 6, 0.8 );

    Sequence fwd1 = Sequence( "fwd1", "ACTGGATGGA", "JJJJJJJJJJ" );
    Sequence rev1 = Sequence( "rev1",      "ATGGAATCCC", "JJJJJJJJJJ" ).ReverseComplement();

    res = merger.Merge( merged, fwd1, rev1 );
    REQUIRE( res == false );
  }

  SECTION( "min identity requirement" ) {
    PairedEndMerger merger( 5, 1.0 );

    Sequence fwd1 = Sequence( "fwd1", "ACTGGATGGA", "JJJJJJJJJJ" );
    Sequence rev1 = Sequence( "rev1",     "GATAGAATCCC", "JJJJJJJJJJ" ).ReverseComplement();

    res = merger.Merge( merged, fwd1, rev1 );
    REQUIRE( res == false );
  }

  SECTION( "posterior Q calculation" ) {
    PairedEndMerger merger( 3, 1.0 );
    Sequence fwd1 = Sequence( "fwd1", "ATTGACCGT",
                                      "1>AA1@FFF" );
    Sequence rev1 = Sequence( "rev1",     "ACCGTGAATC",
                                          "?AAAAFFFFF" ).ReverseComplement();

    res = merger.Merge( merged, fwd1, rev1 );
    REQUIRE( res == true );
    REQUIRE( merged.sequence == "ATTGACCGTGAATC" );
    REQUIRE( merged.quality  == "1>AAJJJJJFFFFF" );
  }
}
