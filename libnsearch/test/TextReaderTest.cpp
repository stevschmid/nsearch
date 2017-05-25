#include <catch.hpp>

#include <nsearch/TextReader.h>

#include <sstream>

TEST_CASE( "TextReader" )  {

#if defined(__APPLE__) || defined(__unix__)
  SECTION( "File" ) {
    const char filename[] = "/tmp/textreadertest.tmp";
    std::ofstream file( filename );
    file << "Hello" << std::endl << std::endl << "Happy " << std::endl << "World";
    file.close();

    std::string line;
    TextFileReader reader( filename );

    reader >> line;
    REQUIRE( line == "Hello" );
    reader >> line;
    REQUIRE( line == "Happy" );
    reader >> line;
    REQUIRE( line == "World" );

    REQUIRE( reader.EndOfFile() == true );

    std::remove( filename );
  }
#endif

  SECTION( "Stream" ) {
    std::istringstream iss( "Hello\nWhat \n\nIs up\n" );

    std::string line;
    TextStreamReader reader( iss );

    reader >> line;
    REQUIRE( line == "Hello" );
    reader >> line;
    REQUIRE( line == "What" );
    reader >> line;
    REQUIRE( line == "Is up" );

    REQUIRE( reader.EndOfFile() == true );
  }
}
