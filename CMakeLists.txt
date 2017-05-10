cmake_minimum_required(VERSION 3.6)
project(lather)

set(CMAKE_CXX_STANDARD 14)
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
        display.cpp)


link_libraries(gsl gslcblas Magick++ sfml-graphics sfml-window sfml-system)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -fno-omit-frame-pointer -march=native -O3 -static -flto -fno-math-errno")

set(CMAKE_BUILD_TYPE Release)

add_executable(display ${SOURCE_FILES})