# This is the CMake file for the src directory of NCEPLIBS-ip.
#
# Mark Potts, Kyle Gerheiser, Ed Hartnett, Eric Engle

# These are the source code filees.
set(fortran_src gdswzd_c.F90 gdswzd_mod.F90 ip2lib_4.h 
ip2lib_d.h ipolates.F90 ipolatev.F90 ipxetas.F90 ipxwafs.F90
ipxwafs2.F90 ipxwafs3.F90 movect.F90 bilinear_interp_mod.F90
bicubic_interp_mod.F90 neighbor_interp_mod.F90 budget_interp_mod.F90
spectral_interp_mod.F90 neighbor_budget_interp_mod.F90 ip_mod.F90
ip_grid_descriptor_mod.F90 ip_grid_mod.F90
ip_station_points_grid_mod.F90 ip_gaussian_grid_mod.F90
ip_equid_cylind_grid_mod.F90 ip_lambert_conf_grid_mod.F90
ip_mercator_grid_mod.F90 ip_polar_stereo_grid_mod.F90
ip_rot_equid_cylind_egrid_mod.F90 ip_rot_equid_cylind_grid_mod.F90
constants_mod.F90 ip_grids_mod.F90 ip_grid_factory_mod.F90
ip_interpolators_mod.F90 earth_radius_mod.F90 polfix_mod.F90)

# Set compiler flags.
if(CMAKE_Fortran_COMPILER_ID MATCHES "^(Intel)$")
  set(CMAKE_Fortran_FLAGS "-g -traceback -warn all -auto -convert big_endian -assume byterecl -fp-model strict -fpp ${CMAKE_Fortran_FLAGS}")
  set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -check all -warn all")
  set(fortran_d_flags "-r8")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "^(GNU)$")
  set(CMAKE_Fortran_FLAGS "-g  -fbacktrace -fconvert=big-endian -cpp ${CMAKE_Fortran_FLAGS}")
  set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -ggdb -Wall -fcheck=all")
  set(fortran_d_flags "-fdefault-real-8")
endif()

# We build a version of the library with 4-byte reals (_4), and one
# with 8-byte reals (_d).
set(kinds "4" "d")
foreach(kind ${kinds})
  set(lib_name ${PROJECT_NAME}_${kind})
  set(module_dir "${CMAKE_CURRENT_BINARY_DIR}/include_${kind}")
  configure_file(ip2lib_${kind}.h ${CMAKE_CURRENT_BINARY_DIR}/include_${kind}/ip2lib.h COPYONLY)
  configure_file(iplib_${kind}.h ${CMAKE_CURRENT_BINARY_DIR}/include_${kind}/iplib.h COPYONLY)

  # Compiled with preprocessor definition LSIZE=D, not d
  string(TOUPPER ${kind} kind_definition)

  set(BUILD_FLAGS "${fortran_${kind}_flags}")

  add_library(${lib_name}_objlib OBJECT ${fortran_src})
  if(BUILD_SHARED_LIBS)
    set_property(TARGET ${lib_name}_objlib PROPERTY POSITION_INDEPENDENT_CODE 1)
  endif()
  add_library(${PROJECT_NAME}::${lib_name} ALIAS ${lib_name}_objlib)

  target_compile_definitions(${lib_name}_objlib PRIVATE "LSIZE=${kind_definition}")
  set_target_properties(${lib_name}_objlib PROPERTIES COMPILE_FLAGS "${BUILD_FLAGS}")
  set_target_properties(${lib_name}_objlib PROPERTIES Fortran_MODULE_DIRECTORY "${module_dir}")
  target_include_directories(${lib_name}_objlib PUBLIC $<BUILD_INTERFACE:${module_dir}> $<INSTALL_INTERFACE:include_${kind}>)

  # Create static object library from object library.
  if(BUILD_STATIC_LIBS)
    add_library(${lib_name}_static STATIC $<TARGET_OBJECTS:${lib_name}_objlib>)
    set_target_properties(${lib_name}_static PROPERTIES OUTPUT_NAME ${lib_name})
    set_target_properties(${lib_name}_static PROPERTIES EXPORT_NAME ${lib_name})
    set_target_properties(${lib_name}_static PROPERTIES EXPORT_NAME ${lib_name}_static)
    target_link_libraries(${lib_name}_static PUBLIC sp::sp_${kind})
    if(OpenMP_Fortran_FOUND)
      target_link_libraries(${lib_name}_static PUBLIC OpenMP::OpenMP_Fortran)
    endif()
    list(APPEND LIB_TARGETS ${lib_name}_static)
  endif()

  # Create shared object library from object library.
  if(BUILD_SHARED_LIBS)
    add_library(${lib_name}_shared SHARED $<TARGET_OBJECTS:${lib_name}_objlib>)
    set_target_properties(${lib_name}_shared PROPERTIES OUTPUT_NAME ${lib_name})
    set_target_properties(${lib_name}_shared PROPERTIES EXPORT_NAME ${lib_name})
    set_target_properties(${lib_name}_shared PROPERTIES SOVERSION 0)
    target_link_libraries(${lib_name}_shared PUBLIC sp::sp_${kind})
    if(OpenMP_Fortran_FOUND)
      target_link_libraries(${lib_name}_shared PUBLIC OpenMP::OpenMP_Fortran)
    endif()
    list(APPEND LIB_TARGETS ${lib_name}_shared)
  endif()
  
  install(DIRECTORY ${module_dir} DESTINATION ${CMAKE_INSTALL_PREFIX})
  install(FILES ip2lib_${kind}.h DESTINATION ${CMAKE_INSTALL_PREFIX}/include_${kind} RENAME ip2lib.h)
  install(FILES iplib_${kind}.h DESTINATION ${CMAKE_INSTALL_PREFIX}/include_${kind} RENAME iplib.h)
endforeach()

# Install the library or libraries.
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
