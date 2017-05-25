#include <catch.hpp>

#include <nsearch/FASTQ/Reader.h>
#include <nsearch/FASTQ/Writer.h>
#include <nsearch/FASTQ/QScore.h>

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

  SECTION( "Q Score" ) {
    const FASTQ::QScore& qscore = FASTQ::QScore::Instance();

    REQUIRE( qscore.ScoreToProbability( 0 ) == 1.0 );
    REQUIRE( qscore.ScoreToProbability( 10 ) == 0.1 );
    REQUIRE( Approx( qscore.ScoreToProbability( 37 ) ).epsilon( 0.00001 ) == 0.00020 );

    REQUIRE( qscore.ProbabilityToScore( 1.0 ) == 0 );
    REQUIRE( qscore.ProbabilityToScore( 0.12589 ) == 9 );
    REQUIRE( qscore.ProbabilityToScore( 0.00008 ) == 41 );

    REQUIRE( qscore.CalculatePosteriorScoreForMatch( 39, 2 ) == 41 );
    REQUIRE( qscore.CalculatePosteriorScoreForMismatch( 39, 2 ) == 38 );

    REQUIRE( qscore.CalculatePosteriorScoreForMatch( 2, 39 ) == 41 );
    REQUIRE( qscore.CalculatePosteriorScoreForMismatch( 2, 39 ) == 38 );

    REQUIRE( qscore.CalculatePosteriorScoreForMatch( 20, 20 ) == std::min( 45, FASTQ::Q_MAX_SCORE ) ); // would be45, but we defined 41 as max QSCORE
    REQUIRE( qscore.CalculatePosteriorScoreForMismatch( 20, 20 ) == 3 );

    REQUIRE( qscore.CalculatePosteriorScoreForMatch( 10, 20 ) == 34 );
    REQUIRE( qscore.CalculatePosteriorScoreForMismatch( 10, 20 ) == 11 ); // 10.508, rounded up

    REQUIRE( qscore.CalculatePosteriorScoreForMatch( 20, 10 ) == 34 );
    REQUIRE( qscore.CalculatePosteriorScoreForMismatch( 20, 10 ) == 11 ); // 10.508, rounded up
  }
}
