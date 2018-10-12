# Specify the minimum version for CMake
cmake_minimum_required(VERSION 2.8)

# Project's name
project(lidCavity)

# Set the output folder where your program will be created
set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})

set_directory_properties(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES ${CMAKE_SOURCE_DIR}/build)

if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
  set (CMAKE_C_FLAGS "-march=native -O3 -std=c11 -Wall")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  set (CMAKE_C_FLAGS "-march=native -flto -O3 -std=gnu11 -Wall")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  set (CMAKE_C_FLAGS "-xHost -O3 -std=c11 -Wall")
endif()

option (USE_OpenMP "Use OpenMP" ON)
if(USE_OpenMP)
  find_package(OpenMP)
  if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  endif()
endif()

add_executable(lidCavity src/lidCavity.c)
target_link_libraries(lidCavity m)
