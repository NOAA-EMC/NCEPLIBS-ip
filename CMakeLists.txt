# This is the main CMake file for NCEPLIBS-ip.
#
# Mark Potts, Kyle Gerheiser, Eric Engle, Ed Hartnett
cmake_minimum_required(VERSION 3.15)

# Get the version from the VERSION file.
file(STRINGS "VERSION" pVersion)

# Set up the cmake project.
project(ip VERSION ${pVersion} LANGUAGES C Fortran)

# Load GNU standard install paths.
include(GNUInstallDirs)

# Handle user options.
option(ENABLE_DOCS "Enable generation of Doxygen-based documentation" OFF)
option(OPENMP "Use OpenMP threading" OFF)
option(BUILD_SHARED_LIBS "Build shared libraries" OFF)
option(BUILD_4 "Build the 4-byte real version of the library, libip_4.{a,so}" ON)
option(BUILD_D "Build the 8-byte real version of the library, libip_d.{a,so}" ON)
option(BUILD_8 "Build the 8-byte integer & real version of the library, libip_8.{a,so}" OFF)
option(BUILD_DEPRECATED "Build deprecated spectral processing functions" OFF)
option(TEST_TIME_LIMIT "Set timeout for tests" OFF)

# Figure whether user wants a _4, a _d, and/or _8.
if(BUILD_4)
  set(kinds "4")
endif()
if(BUILD_D)
  set(kinds ${kinds} "d")
endif()
if(BUILD_8)
  set(kinds ${kinds} "8")
endif()
if(NOT BUILD_4 AND NOT BUILD_D AND NOT BUILD_8)
  message(FATAL_ERROR "At least one of BUILD_4, BUILD_D, and BUILD_8 must be turned on")
endif()
message(STATUS "Library kinds that will be build: ${kinds}")

# Set the build type.
if(NOT CMAKE_BUILD_TYPE MATCHES "^(Debug|Release|RelWithDebInfo|MinSizeRel)$")
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

# Handle OpenMP if desired.
if(OPENMP)
  find_package(OpenMP REQUIRED COMPONENTS Fortran)
endif()

#find_package(LAPACK REQUIRED)

# Set compiler flags.
if(CMAKE_Fortran_COMPILER_ID MATCHES "^(Intel|IntelLLVM)$")
  set(CMAKE_Fortran_FLAGS "-g -traceback -assume byterecl -fp-model strict -fpp -auto ${CMAKE_Fortran_FLAGS}")
  set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -check all -warn all")
  if(CMAKE_Fortran_COMPILER_ID MATCHES "^(IntelLLVM)$")
    # Avoid Intel OneAPI 2023.2.1 bug
    set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -check nouninit ")
  endif()
  set(fortran_d_flags "-r8")
  set(fortran_8_flags "-i8 -r8")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "^(GNU)$")
  set(CMAKE_Fortran_FLAGS "-g -fbacktrace -cpp ${CMAKE_Fortran_FLAGS}")
  set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -ggdb -Wall -Wno-unused-dummy-argument -Wsurprising -Wextra -fcheck=all")
  set(fortran_d_flags "-fdefault-real-8")
  set(fortran_8_flags "-fdefault-integer-8 -fdefault-real-8")
endif()

# This is the source code directiroy.
add_subdirectory(src)

# Build tests.
include(CTest)
if(BUILD_TESTING)
  add_subdirectory(tests)
endif()

# If doxygen documentation we enabled, build it.
if(ENABLE_DOCS)
  find_package(Doxygen REQUIRED)
  set(abs_top_srcdir "${CMAKE_SOURCE_DIR}")
  add_subdirectory(docs)  
endif()
