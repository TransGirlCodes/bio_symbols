cmake_minimum_required(VERSION 3.7)
project(bio_symbols)

set(CMAKE_CXX_STANDARD 14)

set(SOURCE_FILES library.cpp library.h)
add_library(bio_symbols ${SOURCE_FILES})