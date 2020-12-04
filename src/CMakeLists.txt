# This is the CMake file for the src directory of NCEPLIBS-ip.
#
# Mark Potts, Kyle Gerheiser

# These are the source files.
set(fortran_src gausslat.f90 gdswzd00_mod.f90 gdswzd01_mod.f90
    gdswzd03_mod.f90 gdswzd04_mod.f90 gdswzd05_mod.f90 gdswzd_mod.f90
    gdswzdcb_mod.f90 gdswzdcd_mod.f90 ijkgds0.f90 ijkgds1.f90
    ipolates.f90 ipolatev.f90 ipxwafs.f90 ipxwafs2.f90 ipxwafs3.f90
    makgds.f90 movect.f90 polates0.f90 polates1.f90 polates2.f90
    polates3.f90 polates4.f90 polates6.f90 polatev0.f90 polatev1.f90
    polatev2.f90 polatev3.f90 polatev4.f90 polatev6.f90 polfixs.f90
    polfixv.f90 gdswzd_c.F90)

if(CMAKE_Fortran_COMPILER_ID MATCHES "^(Intel)$")
  set(CMAKE_Fortran_FLAGS
    "-g -auto -convert big_endian -assume byterecl -fp-model strict -fpp ${CMAKE_Fortran_FLAGS}"
    )
  set(fortran_d_flags "-i4 -r8")
  set(fortran_8_flags "-i8 -r8")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "^(GNU)$")
  set(CMAKE_Fortran_FLAGS "-g -fconvert=big-endian -cpp ${CMAKE_Fortran_FLAGS}")
  set(CMAKE_Fortran_FLAGS_DEBUG "-ggdb -Wall")
  set(fortran_d_flags "-fdefault-real-8")
  set(fortran_8_flags "-fdefault-integer-8 -fdefault-real-8")
endif()

# We build 3 versions of the library and include files.
set(kinds "4" "8" "d")
foreach(kind ${kinds})
  set(lib_name ${PROJECT_NAME}_${kind})
  set(module_dir "${CMAKE_CURRENT_BINARY_DIR}/include_${kind}")
  configure_file(iplib_${kind}.h ${module_dir}/iplib.h COPYONLY)

  # Compiled with preprocessor definition LSIZE=D, not d
  string(TOUPPER ${kind} kind_definition)

  set(BUILD_FLAGS "${fortran_${kind}_flags}")

  add_library(${lib_name} STATIC ${fortran_src})
  add_library(${PROJECT_NAME}::${lib_name} ALIAS ${lib_name})

  if(OpenMP_Fortran_FOUND)
    target_link_libraries(${lib_name} PRIVATE OpenMP::OpenMP_Fortran)
  endif()

  target_link_libraries(${lib_name} PRIVATE sp::sp_${kind})
  
  target_compile_definitions(${lib_name} PRIVATE "LSIZE=${kind_definition}")
  set_target_properties(${lib_name} PROPERTIES COMPILE_FLAGS "${BUILD_FLAGS}")
  set_target_properties(${lib_name} PROPERTIES Fortran_MODULE_DIRECTORY "${module_dir}")

  target_include_directories(${lib_name}
    PUBLIC $<BUILD_INTERFACE:${module_dir}>
    $<INSTALL_INTERFACE:include_${kind}>)

  list(APPEND LIB_TARGETS ${lib_name})
  
  install(DIRECTORY ${module_dir} DESTINATION ${CMAKE_INSTALL_PREFIX})
endforeach()

install(
  TARGETS ${LIB_TARGETS}
  EXPORT ${PROJECT_NAME}Exports
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})

### Package config
include(CMakePackageConfigHelpers)
set(CONFIG_INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME})

export(EXPORT ${PROJECT_NAME}Exports
  NAMESPACE ${PROJECT_NAME}::
  FILE ${PROJECT_NAME}-targets.cmake)

configure_package_config_file(
  ${CMAKE_SOURCE_DIR}/cmake/PackageConfig.cmake.in ${CMAKE_BINARY_DIR}/${PROJECT_NAME}-config.cmake
  INSTALL_DESTINATION ${CONFIG_INSTALL_DESTINATION})
install(FILES ${CMAKE_BINARY_DIR}/${PROJECT_NAME}-config.cmake
  DESTINATION ${CONFIG_INSTALL_DESTINATION})

write_basic_package_version_file(
  ${CMAKE_BINARY_DIR}/${PROJECT_NAME}-config-version.cmake
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY AnyNewerVersion)
install(FILES ${CMAKE_BINARY_DIR}/${PROJECT_NAME}-config-version.cmake
  DESTINATION ${CONFIG_INSTALL_DESTINATION})

install(EXPORT ${PROJECT_NAME}Exports
  NAMESPACE ${PROJECT_NAME}::
  FILE ${PROJECT_NAME}-targets.cmake
  DESTINATION ${CONFIG_INSTALL_DESTINATION})

