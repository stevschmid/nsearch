# Define library. Only source files here!
project(libnsearch VERSION 0.1 LANGUAGES CXX)

add_library(libnsearch
  src/Sequence.cpp
  src/TextReader.cpp
  src/PairedEnd/Merger.cpp
  )

target_include_directories(libnsearch PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
  PRIVATE src)

target_link_libraries(libnsearch)

# Tests
add_subdirectory(test)
