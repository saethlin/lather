cmake_minimum_required(VERSION 3.6)
project(lather)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_COMPILER /usr/bin/clang++)

set(SOURCE_FILES
        inih/cpp/INIReader.h inih/cpp/INIReader.cpp
        inih/ini.h inih/ini.c
        fitrv.hpp
        simulation.hpp simulation.cpp
        spot.hpp spot.cpp
        star.hpp star.cpp
        profile.hpp profile.cpp
        point.hpp point.cpp
        fitrv.cpp
        boundingshape.cpp boundingshape.hpp
        lather.cpp)

link_libraries(gsl gslcblas)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fno-omit-frame-pointer -march=native -Ofast -Wall -g3")

add_executable(lather ${SOURCE_FILES})