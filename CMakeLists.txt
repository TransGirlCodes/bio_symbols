cmake_minimum_required(VERSION 3.7)

project(biosymbols VERSION 0.1 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 14)
add_compile_options(-Wa -aslh)

add_library(biosymbols INTERFACE)
target_include_directories(biosymbols INTERFACE include/)

add_executable(mytest test/basic_checks.cpp)
target_link_libraries(mytest biosymbols)