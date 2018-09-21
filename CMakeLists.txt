cmake_minimum_required(VERSION 3.12)
project(DQDsystem)

set(CMAKE_CXX_STANDARD 14)

include_directories(DQD)

add_executable(DQDsystem
        DQD/HilbertSpace.cpp
        DQD/HilbertSpace.h
        DQD/MathUtilities.cpp
        DQD/MathUtilities.h
        Test/main.cpp
        DQD/HeitlerLondon.cpp DQD/HeitlerLondon.h)
