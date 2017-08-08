# Define library. Only source files here!
project(libnsearch VERSION 0.1 LANGUAGES CXX)

add_library(libnsearch
  src/TextReader.cpp
  )

target_include_directories(libnsearch PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
  PRIVATE src)

target_link_libraries(libnsearch)

# Zlib
find_package(ZLIB)
if(ZLIB_FOUND)
  add_definitions(-DUSE_ZLIB=1)
  include_directories(${ZLIB_INCLUDE_DIRS})
  target_link_libraries(libnsearch ${ZLIB_LIBRARIES})
endif(ZLIB_FOUND)

# Tests
add_subdirectory(test)
