cmake_minimum_required(VERSION 3.2)
project(auditoryfilterbanks)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")
find_package(Eigen3 REQUIRED)
include_directories(EIGEN3_INCLUDE_DIR)

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

set(SOURCE_FILES
    src/Biquad.cpp
    src/Biquad.h
    src/CochlearFilterbank.cpp
    src/CochlearFilterbank.h
    src/ModulationFilterBank.cpp
    src/ModulationFilterBank.h
    src/Test.cpp)

add_executable(test_auditoryfilterbanks ${SOURCE_FILES})
target_link_libraries(test_auditoryfilterbanks Threads::Threads)

set(CMAKE_BUILD_TYPE Release)
