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
option(ENABLE_DOCS "Enable generation of doxygen-based documentation." OFF)
option(OPENMP "use OpenMP threading" OFF)
option(BUILD_SHARED_LIBS "Build shared libraries" ON)
option(BUILD_STATIC_LIBS "Build static libraries" ON)
option(BUILD_4 "Build the 4-byte real version of the library, libip_4.a" ON)
option(BUILD_D "Build the 8-byte real version of the library, libip_d.a" ON)

# Figure whether user wants a _4, a _d, or both libraries.
if(BUILD_4 AND BUILD_D)
  set(kinds "4" "d")
elseif(BUILD_4 AND NOT BUILD_D)
  set(kinds "4")
elseif(BUILD_D AND NOT BUILD_4)
  set(kinds "d")
else()
  message(FATAL_ERROR "At least one of BUILD_4 or BUILD_D must be turned on")
endif()

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

# We need the NCEPLIBS-sp library.
find_package(sp 2.3.0 REQUIRED)

# If doxygen documentation we enabled, build it.
if(ENABLE_DOCS)
  find_package(Doxygen REQUIRED)
  set(abs_top_srcdir "${CMAKE_SOURCE_DIR}")
  add_subdirectory(docs)  
endif()

# This is the source code directiroy.
add_subdirectory(src)

# Build tests.
include(CTest)
if(BUILD_TESTING)
    add_subdirectory(tests)
endif()
