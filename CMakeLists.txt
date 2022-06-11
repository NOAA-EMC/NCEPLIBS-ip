# This is the main CMake file for NCEPLIBS-ip.
#
# Mark Potts, Kyle Gerheiser
cmake_minimum_required(VERSION 3.15)

# Get the version from the VERSION file.
file(STRINGS "VERSION" pVersion)

project(
  ip
  VERSION ${pVersion}
  LANGUAGES C Fortran)
  
include(GNUInstallDirs)

# Handle build options.
option(ENABLE_DOCS "Enable generation of doxygen-based documentation." OFF)
option(OPENMP "use OpenMP threading" OFF)

if(NOT CMAKE_BUILD_TYPE MATCHES "^(Debug|Release|RelWithDebInfo|MinSizeRel)$")
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

# We only have valid options for the 4,8,d builds with these compilers.
if(NOT CMAKE_Fortran_COMPILER_ID MATCHES "^(Intel|GNU)$")
  message(WARNING "Compiler not officially supported: ${CMAKE_Fortran_COMPILER_ID}")
endif()

# If the user wants to use OpenMP, find it.
if(OPENMP)
  find_package(OpenMP REQUIRED COMPONENTS Fortran)
endif()

find_package(sp 2.3.0 REQUIRED)

# If doxygen documentation we enabled, build it. This must come before
# adding the source code directories; the main documentation build
# must happen before any of the utility document builds.
if(ENABLE_DOCS)
  find_package(Doxygen REQUIRED)
  set(abs_top_srcdir "${CMAKE_SOURCE_DIR}")
  add_subdirectory(docs)  
endif()

# Build the library.
add_subdirectory(src)

# Build and run tests.
include(CTest)
if(BUILD_TESTING)
    add_subdirectory(tests)
endif()
