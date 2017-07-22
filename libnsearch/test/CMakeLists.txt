enable_testing()

add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND})

set(CATCH_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/vendor/catch)
add_library(Catch INTERFACE)
target_include_directories(Catch INTERFACE ${CATCH_INCLUDE_DIR})

add_executable(testnsearch EXCLUDE_FROM_ALL
  Test.cpp
  SequenceTest.cpp
  TextReaderTest.cpp
  FASTQTest.cpp
  FASTATest.cpp
  PairedEndTest.cpp
  UtilsTest.cpp
  Alignment/RangesTest.cpp
  Alignment/HitTrackerTest.cpp
  Alignment/ExtendAlignTest.cpp
  Alignment/BandedAlignTest.cpp
  Database/KmersTest.cpp
  Database/HSPTest.cpp
  )

add_dependencies(check testnsearch)

target_link_libraries(testnsearch
  libnsearch
  Catch)

add_test(NAME testnsearch COMMAND testnsearch)

# Parse tests
include(${CATCH_INCLUDE_DIR}/ParseAndAddCatchTests.cmake)
ParseAndAddCatchTests(testnsearch)
