add_executable(nsearch
  src/Main.cpp)

set(DOCOPT_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/vendor/docopt)
add_library(docopt INTERFACE)
target_include_directories(docopt INTERFACE ${DOCOPT_INCLUDE_DIR})
target_compile_definitions(docopt INTERFACE DOCOPT_HEADER_ONLY=1)

target_link_libraries(nsearch
  libnsearch
  docopt)
