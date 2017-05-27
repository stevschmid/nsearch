#include <catch.hpp>

#include <nsearch/Aligner.h>

TEST_CASE( "Aligner" )  {
  Aligner aligner( 2, -3, 5, 2 );

  SECTION( "Global" ) {
    Sequence seq1( "GACTTAC");
    Sequence seq2( "CGTGAATTCAT" );

    SECTION( "Default" ) {
      GlobalAlignment aln = aligner.GlobalAlign( seq1, seq2 );
      REQUIRE( aln.score == -14 );
      REQUIRE( CigarAsString( aln.cigar ) == "3D5M1D2M" );
    }

    SECTION( "Symmetry" ) {
      REQUIRE( aligner.GlobalAlign( seq1, seq2 ).score == aligner.GlobalAlign( seq2, seq1 ).score );
    }
  }

  SECTION( "Local" ) {
    Sequence seq1( "TACGGGCCCGCTAC" );
    Sequence seq2( "TAGCCCTATCGGTCA" );

    SECTION( "Gaps" ) {
      LocalAlignment aln1 = Aligner( 5, -4, 4, 1 ).LocalAlign( seq1, seq2 );
      REQUIRE( aln1.score == 27 );


      LocalAlignment aln2 = Aligner( 5, -4, 10, 1 ).LocalAlign( seq1, seq2 );
      REQUIRE( aln2.score == 20 );
    }

    SECTION( "Ambiguous Nucleotides" ) {
      Sequence seq1 = "ATGCTGGTACCTGGGAT";
      Sequence seq2 =     "TGGYAWNNV"; // W matches A or T, here C (mismatch)
      LocalAlignment aln = aligner.LocalAlign( seq1, seq2 );
      REQUIRE( aln.score == 8 * 2 - 1 * 3 );
    }

    SECTION( "Symmetry" ) {
      Aligner aligner;
      REQUIRE( aligner.LocalAlign( seq1, seq2 ).score == aligner.LocalAlign( seq2, seq1 ).score );
    }
  }

  /* SECTION( "Reader" ) { */
  /*   std::string content = ">Seq1\n" */
  /*     "TGGCG\n" */
  /*     "ATTGG\n" */
  /*     "\n" */
  /*     ">Seq2\n" */
  /*     "TTTTT\n" */
  /*     "CAGTC\n" */
  /*     ">Seq3\n" */
  /*     "actgc\n"; */

  /*   std::istringstream iss( content ); */

  /*   FASTA::Reader reader( iss ); */
  /*   Sequence sequence; */

  /*   reader >> sequence; */
  /*   REQUIRE( sequence.identifier == "Seq1" ); */
  /*   REQUIRE( sequence.sequence == "TGGCGATTGG" ); */

  /*   reader >> sequence; */
  /*   REQUIRE( sequence.identifier == "Seq2" ); */
  /*   REQUIRE( sequence.sequence == "TTTTTCAGTC" ); */

  /*   reader >> sequence; */
  /*   REQUIRE( sequence.identifier == "Seq3" ); */
  /*   REQUIRE( sequence.sequence == "ACTGC" ); */

  /*   REQUIRE( reader.EndOfFile() == true ); */
  /* } */
}
