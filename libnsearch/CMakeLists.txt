# Define library. Only source files here!
project(libnsearch VERSION 0.1 LANGUAGES CXX)

add_library(libnsearch
  src/Sequence.cpp)

target_include_directories(libnsearch PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
  PRIVATE src)

# Test
set(CATCH_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/vendor/catch)
add_library(Catch INTERFACE)
target_include_directories(Catch INTERFACE ${CATCH_INCLUDE_DIR})

add_executable(testnsearch
  test/Test.cpp
  test/SequenceTest.cpp)

target_link_libraries(testnsearch
  libnsearch
  Catch)

add_test(NAME testnsearch COMMAND testnsearch)
