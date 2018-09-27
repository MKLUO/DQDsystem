cmake_minimum_required(VERSION 3.12)
project(DQDSystem)

set(CMAKE_CXX_STANDARD 14)

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
target_link_libraries(DQDSystem Fourier DQD Plot)

#------ArrayFire & CUDA
#
#find_package(ArrayFire REQUIRED)
#include_directories(${ArrayFire_INCLUDE_DIRS})
#
## ArrayFire CUDA backend
#find_package(CUDA)
#IF(${ArrayFire_CUDA_FOUND} AND ${CUDA_FOUND})
#    # We need to find CUDA and NVVM as transitive linking is disabled on some OSes
#    find_package(CUDA REQUIRED)
#    #find_package(NVVM REQUIRED)
#    set(NVVM_LIB
#            "C:/Program Files/NVIDIA GPU Computing Toolkit/CUDA/v10.0/nvvm/lib/x64/nvvm.lib")
#    #set(ArrayFire_CUDA_LIBRARIES "C:/Program Files/ArrayFire/v3/lib/afcuda.lib")
#    message(STATUS ${CUDA_TOOLKIT_ROOT_DIR})
#    message(STATUS "ArrayFire CUDA found. Enabling CUDA benchmark")
#    message(STATUS ${ArrayFire_CUDA_LIBRARIES})
#    message(STATUS ${CMAKE_THREAD_LIBS_INIT})
#    message(STATUS ${CUDA_LIBRARIES})
#    message(STATUS ${NVVM_LIB})
#    target_link_libraries(Fourier ${ArrayFire_CUDA_LIBRARIES}
#            ${CMAKE_THREAD_LIBS_INIT} ${CUDA_LIBRARIES} ${NVVM_LIB})
#ENDIF()