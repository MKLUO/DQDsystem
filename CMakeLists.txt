cmake_minimum_required(VERSION 3.12)
project(DQDSystem)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories(DQD
        Utilities)

add_library(Fourier
        Utilities/fourier.hpp
        Utilities/fourier.cpp)

add_library(DQD
        DQD/hilbert_space.cpp
        DQD/hilbert_space.hpp

        DQD/math_utilities.cpp
        DQD/math_utilities.hpp

        DQD/hl_hm.cpp
        DQD/hl_hm.hpp)

add_library(Plot
        Utilities/plot.hpp
        Utilities/plot.cpp)

add_executable(DQDSystem
        Test/test.cpp)

target_link_libraries(DQD Fourier)

#------ArrayFire & CUDA

find_package(ArrayFire)
IF(${ArrayFire_FOUND})
    include_directories(${ArrayFire_INCLUDE_DIRS})
ELSE()
    message(FATAL_ERROR "ArrayFire not found. Please install ArrayFire.")
ENDIF()

find_package(CUDA)
# ArrayFire CUDA/CPU backend
# set(CUDA_FOUND FALSE)
IF(${CUDA_FOUND})
    message(STATUS "ArrayFire CUDA found. Enabling CUDA/CPU benchmark")

    target_link_libraries(Fourier ${ArrayFire_LIBRARIES}
            ${CMAKE_THREAD_LIBS_INIT} ${CUDA_LIBRARIES})

    target_link_libraries(DQDSystem Fourier DQD Plot ${ArrayFire_LIBRARIES}
            ${CMAKE_THREAD_LIBS_INIT} ${CUDA_LIBRARIES})

# ArrayFire CPU backend
ELSE()
    message(STATUS "ArrayFire CUDA not found. Using CPU backend.")

    target_link_libraries(Fourier ${ArrayFire_CPU_LIBRARIES}
            ${CMAKE_THREAD_LIBS_INIT})

    target_link_libraries(DQDSystem Fourier DQD Plot ${ArrayFire_CPU_LIBRARIES}
            ${CMAKE_THREAD_LIBS_INIT})

ENDIF()

