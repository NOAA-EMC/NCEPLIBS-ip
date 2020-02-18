cmake_minimum_required(VERSION 3.15)

project(
  ip
  VERSION 3.0.0
  LANGUAGES Fortran)

if(NOT CMAKE_BUILD_TYPE MATCHES "^(Debug|Release|RelWithDebInfo|MinSizeRel)$")
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE
      "Release"
      CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
                                               "MinSizeRel" "RelWithDebInfo")
endif()

if(NOT CMAKE_Fortran_COMPILER_ID MATCHES "^(Intel|GNU|Clang|AppleClang|PGI)$")
  message(
    WARNING "Compiler not officially supported: ${CMAKE_Fortran_COMPILER_ID}")
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES "^(Intel)$")
  set(CMAKE_Fortran_FLAGS
      "-auto -qopenmp -convert big_endian -assume byterecl -fp-model strict -fpp"
  )
  set(CMAKE_Fortran_FLAGS_DEBUG "-g")
  set(fortran_d_flag "-i4 -r8")
  set(fortran_8_flag "-i8 -r8")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$")
  set(CMAKE_Fortran_FLAGS "-fconvert=big-endian -cpp")
  set(CMAKE_Fortran_FLAGS_DEBUG "-g")
  set(fortran_d_flags "-fdefault-real-8")
  set(fortran_8_flags "-fdefault-integer-8 -fdefault-real-8")
endif()

add_subdirectory(src)

set(kinds "4" "8" "d")
foreach(kind ${kinds})
  set(lib_name ${PROJECT_NAME}_${kind})
  set(versioned_lib_name ${PROJECT_NAME}_v${PROJECT_VERSION}_${kind})
  set(module_dir "${CMAKE_CURRENT_BINARY_DIR}/include_${kind}")

  # Compiled with preprocessor definition LSIZE=D, not d
  string(TOUPPER ${kind} kind_definition)

  add_library(${lib_name} STATIC ${fortran_src})
  target_compile_definitions(${lib_name} PRIVATE "LSIZE=${kind_definition}")
  set_target_properties(${lib_name} PROPERTIES COMPILE_FLAGS
                                               "${fortran_${kind}_flags}")
  set_target_properties(${lib_name} PROPERTIES OUTPUT_NAME
                                               ${versioned_lib_name})
  set_target_properties(${lib_name} PROPERTIES Fortran_MODULE_DIRECTORY
                                               "${module_dir}")
  target_include_directories(
    ${lib_name}
    PUBLIC $<BUILD_INTERFACE:${module_dir}>
           $<INSTALL_INTERFACE:${CMAKE_INSTALL_PREFIX}/include_${kind}>)

  install(
    TARGETS ${lib_name}
    EXPORT ${PROJECT_NAME}-config
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib)
  install(DIRECTORY ${module_dir} DESTINATION ${CMAKE_INSTALL_PREFIX})
endforeach()

install(EXPORT ${PROJECT_NAME}-config
        DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/cmake/${PROJECT_NAME})
