cmake_minimum_required(VERSION 3.1.0 FATAL_ERROR)
project(nsearch VERSION 0.1 LANGUAGES CXX)

# Enable C++11
set(CMAKE_CXX_STANDARD 11)

# Targets
add_subdirectory(libnsearch)
add_subdirectory(app)
add_subdirectory(vendor/docopt)
