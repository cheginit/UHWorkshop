cmake_minimum_required(VERSION 3.9)

project(lidCavity
    VERSION 0.3
    LANGUAGES C
)

add_executable(lidCavity
    src/lidCavity.c
)

set_directory_properties(
    PROPERTIES
    ADDITIONAL_MAKE_CLEAN_FILES ${CMAKE_SOURCE_DIR}/build
)

set_target_properties(lidCavity
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin
)

target_compile_options(lidCavity
  PUBLIC
  -march=native -O3 -std=c11 -Wall
)

if("${CMAKE_C_COMPILER_ID}" STREQUAL "Intel")
  target_compile_options(lidCavity
    PUBLIC
    -xHost -O3 -std=c11 -Wall
  )
endif()

option (USE_OpenMP "Use OpenMP" ON)
if(USE_OpenMP)
  find_package(OpenMP)
  if(OPENMP_FOUND)
    set(OMP_LIB "OpenMP::OpenMP_C")
  endif()
endif()

target_link_libraries(lidCavity
  PUBLIC
  ${OMP_LIB}
  m
)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
