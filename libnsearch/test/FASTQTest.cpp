#include <catch.hpp>

#include <nsearch/FASTQ/Reader.h>
#include <nsearch/FASTQ/Writer.h>

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

  SECTION( "Writer" ) {
    Sequence seq1( "Seq1", "TAGGC", "JJ:BB" );
    Sequence seq2( "Seq2", "CTAGG", "AA..D" );

    std::ostringstream oss;
    FASTQ::Writer writer( oss );

    writer << seq1;
    writer << seq2;

    const char *expectedOutput = "@Seq1\n"
                                 "TAGGC\n"
                                 "+\n"
                                 "JJ:BB\n"
                                 "@Seq2\n"
                                 "CTAGG\n"
                                 "+\n"
                                 "AA..D\n";
    REQUIRE( oss.str() == expectedOutput );
  }
}
