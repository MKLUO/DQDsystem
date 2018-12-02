cmake_minimum_required(VERSION 3.12)
project(DQDSystem)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_library(Math
        src/Utilities/math_utilities.cpp
        src/Utilities/fourier.cpp)

add_library(DQD
        src/DQD/hilbert_space.cpp
        src/DQD/hl_hm.cpp)

add_library(Plot
        src/Utilities/plot.cpp)

target_link_libraries(DQD Math)
target_link_libraries(Plot Math)

target_include_directories(DQD PUBLIC 
        src/Utilities
        src/DQD)

#------ArrayFire & CUDA

find_package(ArrayFire)
IF(${ArrayFire_FOUND})
    include_directories(${ArrayFire_INCLUDE_DIRS})
ELSE()
    message(FATAL_ERROR "ArrayFire not found. Please install ArrayFire.")
ENDIF()

find_package(CUDA)
# ArrayFire CUDA/CPU backend
IF(${CUDA_FOUND})
    message(STATUS "ArrayFire CUDA found. Enabling CUDA/CPU benchmark")

    target_link_libraries(Math ${ArrayFire_LIBRARIES}
            ${CMAKE_THREAD_LIBS_INIT} ${CUDA_LIBRARIES})

# ArrayFire CPU backend
ELSE()
    message(STATUS "ArrayFire CUDA not found. Using CPU backend.")

    target_link_libraries(Math ${ArrayFire_CPU_LIBRARIES}
            ${CMAKE_THREAD_LIBS_INIT})

ENDIF()

