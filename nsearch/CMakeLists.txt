cmake_minimum_required(VERSION 3.1.0 FATAL_ERROR)

# Setup project
set(APP_NAME nsearch)
set(APP_VERSION 0.1)

project(${APP_NAME} VERSION ${APP_VERSION} LANGUAGES CXX)

# Enable C++11
set(CMAKE_CXX_STANDARD 11)

# We need some defines
add_definitions(-DAPP_NAME="${APP_NAME}")
add_definitions(-DAPP_VERSION=${APP_VERSION})
string(TIMESTAMP BUILD_TIMESTAMP %Y-%m-%d)
add_definitions(-DBUILD_TIMESTAMP="${BUILD_TIMESTAMP}")

# Targets
add_executable(nsearch
  src/Main.cpp
  src/Merge.cpp
  src/Search.cpp
  src/Filter.cpp
  )

add_subdirectory(vendor/docopt)

target_link_libraries(nsearch
  libnsearch
  docopt)
