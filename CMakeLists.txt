cmake_minimum_required(VERSION 3.7)
project(biosymbols)

set(CMAKE_CXX_STANDARD 14)

set(SOURCE_FILES library.cpp library.h)
add_library(biosymbols ${SOURCE_FILES})