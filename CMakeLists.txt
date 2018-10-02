cmake_minimum_required(VERSION 3.12)
project(DQDSystem)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories(DQD
        Utilities)

add_library(Fourier
        Utilities/Fourier.h
        Utilities/Fourier.cpp)

add_library(DQD
        DQD/HilbertSpace.cpp
        DQD/HilbertSpace.h

        DQD/MathUtilities.cpp
        DQD/MathUtilities.h

        DQD/HeitlerLondon.cpp
        DQD/HeitlerLondon.h)

add_library(Plot
        Utilities/Plot.h
        Utilities/Plot.cpp)

add_executable(DQDSystem
        Test/main.cpp)

target_link_libraries(DQD Fourier)

#------ArrayFire & CUDA

find_package(ArrayFire)
IF(${ArrayFire_FOUND})
    include_directories(${ArrayFire_INCLUDE_DIRS})
ELSE()
    message(FATAL_ERROR "ArrayFire not found. Please install ArrayFire.")
ENDIF()

find_package(CUDA)
# ArrayFire CUDA backend
# set(CUDA_FOUND FALSE)
IF(${ArrayFire_CUDA_FOUND} AND ${CUDA_FOUND})
    message(STATUS "ArrayFire CUDA found. Enabling CUDA benchmark")

    target_link_libraries(Fourier ${ArrayFire_CUDA_LIBRARIES}
            ${CMAKE_THREAD_LIBS_INIT} ${CUDA_LIBRARIES})

    target_link_libraries(DQDSystem Fourier DQD Plot ${ArrayFire_CUDA_LIBRARIES}
            ${CMAKE_THREAD_LIBS_INIT} ${CUDA_LIBRARIES})

# ArrayFire CPU backend
ELSE()
    message(STATUS "ArrayFire CUDA not found. Using CPU backend.")

    target_link_libraries(Fourier ${ArrayFire_CPU_LIBRARIES}
            ${CMAKE_THREAD_LIBS_INIT})

    target_link_libraries(DQDSystem Fourier DQD Plot ${ArrayFire_CPU_LIBRARIES}
            ${CMAKE_THREAD_LIBS_INIT})

ENDIF()

