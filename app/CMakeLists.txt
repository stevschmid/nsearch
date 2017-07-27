add_executable(nsearch
  src/Main.cpp
  src/Merge.cpp
  src/Search.cpp
  )

target_link_libraries(nsearch
  libnsearch
  docopt)
