#include <catch.hpp>

#include <nsearch/Utils.h>

TEST_CASE( "Utils" ) {
  SECTION( "UpcaseString" ) {
    std::string str = "AcGt";
    UpcaseString( str );
    REQUIRE( str == "ACGT" );
  }
}
