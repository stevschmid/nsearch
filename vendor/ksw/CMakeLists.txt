project(ksw)
add_library(ksw ksw.c)
target_include_directories(ksw INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
