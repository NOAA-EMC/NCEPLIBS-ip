# This is the main CMake file for NCEPLIBS-ip.
#
# Mark Potts, Kyle Gerheiser
cmake_minimum_required(VERSION 3.15)

# Get the version from the VERSION file.
file(STRINGS "VERSION" pVersion)

project(
  ip
  VERSION ${pVersion}
  LANGUAGES Fortran)

option(OPENMP "use OpenMP threading" OFF)
option(ENABLE_TESTS "Build pfunit tests?" OFF)

if(NOT CMAKE_BUILD_TYPE MATCHES "^(Debug|Release|RelWithDebInfo|MinSizeRel)$")
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE
    "Release"
    CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
    "MinSizeRel" "RelWithDebInfo")
endif()

if(NOT CMAKE_Fortran_COMPILER_ID MATCHES "^(Intel|GNU|Clang|AppleClang)$")
  message(
    WARNING "Compiler not officially supported: ${CMAKE_Fortran_COMPILER_ID}")
endif()

if(OPENMP)
  find_package(OpenMP REQUIRED COMPONENTS Fortran)
endif()

# The NCEPLIBS-sp library is required.
find_package(sp REQUIRED)

add_subdirectory(src)

if (ENABLE_TESTS)
  find_package(PFUNIT REQUIRED)
  enable_testing()
  add_subdirectory(tests)
endif()

# Determine whether or not to generate documentation.
option(ENABLE_DOCS "Enable generation of doxygen-based documentation." OFF)
if(ENABLE_DOCS)
  find_package(Doxygen REQUIRED)
endif()
add_subdirectory(docs)  
