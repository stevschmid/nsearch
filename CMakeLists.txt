cmake_minimum_required(VERSION 3.1.0 FATAL_ERROR)
project(nsearch VERSION 0.1 LANGUAGES CXX)

# Output folders
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Enable C++11
set(CMAKE_CXX_STANDARD 11)

# Targets
enable_testing()
add_subdirectory(libnsearch)
add_subdirectory(app)
add_subdirectory(vendor/docopt)
