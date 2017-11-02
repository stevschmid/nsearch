#include <catch.hpp>

#include <nsearch/CSV/Writer.h>

#include <sstream>

const char CSVOutput[] = R"(QueryId,TargetId,QueryMatchStart,QueryMatchEnd,TargetMatchStart,TargetMatchEnd,QueryMatchSeq,TargetMatchSeq,NumColumns,NumMatches,NumMismatches,NumGaps,Identity,Alignment
"Query,1","Ref,1",1,16,4,22,ATCGTGTACCAGGATG,ATCGTGTCCCACCAGGATG,19,16,0,3,0.842,7=3D9=
"Query,1",Ref2,16,2,3,17,CATCCTGGTACACGA,CATCCTCGTACACGA,15,14,1,0,0.933,6=1X8=
)";

TEST_CASE( "CSV" ) {
  auto entry = std::make_pair(
    Sequence< DNA >( "Query,1", "ATCGTGTACCAGGATG" ),
    HitList< DNA >( {
      { { "Ref,1", "TTTATCGTGTCCCACCAGGATGTTT" }, "3D7=3D9=3D", DNA::Strand::Plus },
      /*
       *
       *   ATCGTGTACCAGGATG (Query +Strand)
       *   CATCCTGGTACACGAT (Query -Strand)
       *   |||||| ||||||||
       * TTCATCCTCGTACACGA- (Database +Strand)
       *
       */
      { { "Ref2", "TTCATCCTCGTACACGA" }, "2D6=1X8=1I", DNA::Strand::Minus },
    } ) );

  std::ostringstream oss;
  CSV::Writer< DNA > writer( oss );

  writer << entry;
  REQUIRE( oss.str() == CSVOutput );
}
