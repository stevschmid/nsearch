#include <catch.hpp>

#include <nsearch/FASTA/Reader.h>
#include <nsearch/FASTA/Writer.h>
#include <nsearch/Alphabet/DNA.h>
#include <nsearch/Alphabet/Protein.h>

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

  SECTION( "Writer" ) {
    Sequence< DNA > seq1( "Seq1", "MTEITAAMVKELRESTGAGMMDCKNALSETNGDFDKAVQLLREKGLGKAAKKADRLAAEGLVSVKVSDDFTIAAMRPSYLSYEDLDMTFVENEYKALVAELEKENEERRRL" );
    Sequence< DNA > seq2( "Seq2", "SATVSEINSETDFVAKNDQFIALTKDTTAHIQSNSLQSVEELHSSTINGVKFEEYLKSQI" );
    Sequence< DNA > seq3( "Seq3", "IALT" );

    std::ostringstream oss;
    FASTA::Writer< DNA > writer( oss );

    writer << seq1 << seq2 << seq3;

    const char* expectedOutput = ">Seq1\n"
                                 "MTEITAAMVKELRESTGAGMMDCKNALSETNGDFDKAVQLLREKGLGKAAKKADRLAAEG\n"
                                 "LVSVKVSDDFTIAAMRPSYLSYEDLDMTFVENEYKALVAELEKENEERRRL\n"
                                 ">Seq2\n"
                                 "SATVSEINSETDFVAKNDQFIALTKDTTAHIQSNSLQSVEELHSSTINGVKFEEYLKSQI\n"
                                 ">Seq3\n"
                                 "IALT\n";
    REQUIRE( oss.str() == expectedOutput );
  }
}
