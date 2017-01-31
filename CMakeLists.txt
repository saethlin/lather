cmake_minimum_required(VERSION 3.6)
project(lather)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_COMPILER /usr/bin/clang++)

set(SOURCE_FILES
        inih/cpp/INIReader.cpp
        inih/cpp/INIReader.h
        inih/ini.c
        inih/ini.h
        fitrv.cpp
        fitrv.hpp
        fitsim.cpp
        fitsim.hpp
        lather.cpp
        profile.cpp
        profile.hpp
        simulation.cpp
        simulation.hpp
        spot.cpp
        spot.hpp
        star.cpp
        star.hpp)

link_libraries(gsl gslcblas)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -march=native -fno-omit-frame-pointer -O3 -funroll-loops")

add_executable(lather ${SOURCE_FILES})