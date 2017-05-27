# Define library. Only source files here!
project(libnsearch VERSION 0.1 LANGUAGES CXX)

add_library(libnsearch
  src/Sequence.cpp
  src/TextReader.cpp
  src/PairedEnd/Merger.cpp
  src/LocalAligner.cpp
  )

target_include_directories(libnsearch PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
  PRIVATE src)

target_link_libraries(libnsearch
  ksw)

# Tests
add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND})

set(CATCH_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/vendor/catch)
add_library(Catch INTERFACE)
target_include_directories(Catch INTERFACE ${CATCH_INCLUDE_DIR})

add_executable(testnsearch EXCLUDE_FROM_ALL
  test/Test.cpp
  test/SequenceTest.cpp
  test/TextReaderTest.cpp
  test/FASTQTest.cpp
  test/FASTATest.cpp
  test/PairedEndTest.cpp
  )

add_dependencies(check testnsearch)

target_link_libraries(testnsearch
  libnsearch
  Catch)

add_test(NAME testnsearch COMMAND testnsearch)

# Parse tests
include(${CATCH_INCLUDE_DIR}/ParseAndAddCatchTests.cmake)
ParseAndAddCatchTests(testnsearch)
