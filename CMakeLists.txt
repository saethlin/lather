cmake_minimum_required(VERSION 3.6)
project(lather)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_COMPILER /usr/bin/clang++)

set(SOURCE_FILES
        inih/cpp/INIReader.h inih/cpp/INIReader.cpp
        inih/ini.h inih/ini.c
        fitrv.hpp fitrv.cpp
        compute_bisector.hpp compute_bisector.cpp
        simulation.hpp simulation.cpp
        spot.hpp spot.cpp
        star.hpp star.cpp
        profile.hpp profile.cpp
        point.hpp point.cpp
        boundingshape.cpp boundingshape.hpp
        lather.cpp
        python_interface.cpp)

include_directories(/home/ben/anaconda3/include/python3.6m)
link_directories(/home/ben/anaconda3/lib)

#link_libraries(gsl gslcblas)
link_libraries(gsl gslcblas python3.6m)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -Wall -Wno-writable-strings -O3")
set(CMAKE_EXE_LINKER_FLAGS -s)

add_library(lather SHARED ${SOURCE_FILES})
#add_executable(lather ${SOURCE_FILES})
