enable_testing()

add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND})

set(CATCH_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/vendor/catch)
add_library(Catch INTERFACE)
target_include_directories(Catch INTERFACE ${CATCH_INCLUDE_DIR})

add_executable(testnsearch EXCLUDE_FROM_ALL
  Alignment/BandedAlignTest.cpp
  Alignment/CigarTest.cpp
  Alignment/ExtendAlignTest.cpp
  Alnout/WriterTest.cpp
  Alphabet/DNATest.cpp
  Alphabet/ProteinTest.cpp
  Database/GlobalSearchTest.cpp
  Database/HSPTest.cpp
  Database/KmersTest.cpp
  DatabaseTest.cpp
  FASTATest.cpp
  FASTQTest.cpp
  PairedEndTest.cpp
  SequenceTest.cpp
  Test.cpp
  TextReaderTest.cpp
  UtilsTest.cpp
  )

add_dependencies(check testnsearch)

target_link_libraries(testnsearch
  libnsearch
  Catch)

add_test(NAME testnsearch COMMAND testnsearch)

# Parse tests
include(${CATCH_INCLUDE_DIR}/ParseAndAddCatchTests.cmake)
ParseAndAddCatchTests(testnsearch)
