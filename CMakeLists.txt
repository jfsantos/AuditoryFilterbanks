cmake_minimum_required(VERSION 2.8)
project(auditoryfilterbanks)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread -std=c++11")

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")
find_package(Eigen3 REQUIRED)
include_directories(EIGEN3_INCLUDE_DIR)

MESSAGE( STATUS "EIGEN3_INCLUDE_PATH: " ${EIGEN3_INCLUDE_PATH} )

set(SOURCE_FILES
    CochlearFilterbank/Biquad.cpp
    CochlearFilterbank/Biquad.h
    CochlearFilterbank/CochlearFilterbank.cpp
    CochlearFilterbank/CochlearFilterbank.h
    CochlearFilterbank/ModulationFilterBank.cpp
    CochlearFilterbank/ModulationFilterBank.h
    CochlearFilterbank/Test.cpp)

add_executable(auditoryfilterbanks ${SOURCE_FILES})
