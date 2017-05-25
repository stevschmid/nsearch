#include <catch.hpp>

#include <nsearch/FASTQ/Reader.h>

#include <sstream>

TEST_CASE( "FASTQ" )  {
  SECTION( "Reader" ) {
    std::string content = "@Seq1\n"
      "TGGCG\n"
      "+\n"
      "JJJJB\n"
      "@Seq2 \n"
      "actgc\n"
      "+\n"
      "JAJI=";

    std::istringstream iss( content );

    FASTQ::Reader reader( iss );
    Sequence sequence;

    reader >> sequence;
    REQUIRE( sequence.identifier == "Seq1" );
    REQUIRE( sequence.sequence == "TGGCG" );
    REQUIRE( sequence.quality == "JJJJB" );

    reader >> sequence;
    REQUIRE( sequence.identifier == "Seq2" );
    REQUIRE( sequence.sequence == "ACTGC" );
    REQUIRE( sequence.quality == "JAJI=" );

    REQUIRE( reader.EndOfFile() == true );
  }
}
