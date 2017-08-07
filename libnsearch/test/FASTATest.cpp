#include <catch.hpp>

#include <nsearch/FASTA/Reader.h>
#include <nsearch/Alphabet/DNA.h>

#include <sstream>

TEST_CASE( "FASTA" ) {
  SECTION( "Reader" ) {
    std::string content = ">Seq1\n"
                          "TGGCG\n"
                          "ATTGG\n"
                          "\n"
                          ">Seq2\n"
                          "TTTTT\n"
                          "CAGTC\n"
                          ">Seq3\n"
                          "actgc\n";

    std::istringstream iss( content );

    FASTA::Reader< DNA > reader( iss );
    Sequence< DNA > sequence;

    reader >> sequence;
    REQUIRE( sequence.identifier == "Seq1" );
    REQUIRE( sequence.sequence == "TGGCGATTGG" );

    reader >> sequence;
    REQUIRE( sequence.identifier == "Seq2" );
    REQUIRE( sequence.sequence == "TTTTTCAGTC" );

    reader >> sequence;
    REQUIRE( sequence.identifier == "Seq3" );
    REQUIRE( sequence.sequence == "ACTGC" );

    REQUIRE( reader.EndOfFile() == true );
  }
}
