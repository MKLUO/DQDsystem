cmake_minimum_required(VERSION 3.12)
project(DQDsystem)

set(CMAKE_CXX_STANDARD 14)

include_directories(DQD)

add_executable(DQDsystem
        DQD/DQD.cpp
        DQD/DQD.h
        DQD/MathUtilities.cpp
        DQD/MathUtilities.h
        Test/main.cpp
        )
