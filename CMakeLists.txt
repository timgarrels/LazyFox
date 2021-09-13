cmake_minimum_required(VERSION 3.16)
project(FOX)

#set(CMAKE_CXX_STANDARD 20 -fopenmp)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++2a -fopenmp")

set(SOURCE_FILES main.cpp tqdm.hpp json.hpp)

add_executable(LazyFox  ${SOURCE_FILES})

set(CMAKE_CXX_FLAGS_RELEASE "-O3")
