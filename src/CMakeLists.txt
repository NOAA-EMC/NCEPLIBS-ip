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
foreach(kind ${kinds})
  set(lib_name ${PROJECT_NAME}_${kind})
  set(module_dir "${CMAKE_CURRENT_BINARY_DIR}/include_${kind}")

  # Create the C header files for _4 or _d.
  configure_file(ip2lib_${kind}.h ${CMAKE_CURRENT_BINARY_DIR}/include_${kind}/ip2lib.h COPYONLY)
  configure_file(iplib_${kind}.h ${CMAKE_CURRENT_BINARY_DIR}/include_${kind}/iplib.h COPYONLY)

  # Compiled with preprocessor definition LSIZE=D, not d
  string(TOUPPER ${kind} kind_definition)

  # Build the library.
  add_library(${lib_name} ${fortran_src})
  add_library(${PROJECT_NAME}::${lib_name} ALIAS ${lib_name})
  
  # For shared, set SOVERSION; increment it if API changes (other than
  # adding new functions).
  if(BUILD_SHARED_LIBS)
    set_target_properties(${lib_name} PROPERTIES SOVERSION 0)
  endif()

  # Set compiler flags.
  target_compile_definitions(${lib_name} PRIVATE "LSIZE=${kind_definition}")
  set(BUILD_FLAGS "${fortran_${kind}_flags}")
  set_target_properties(${lib_name} PROPERTIES COMPILE_FLAGS "${BUILD_FLAGS}")
  set_target_properties(${lib_name} PROPERTIES Fortran_MODULE_DIRECTORY "${module_dir}")
  target_include_directories(${lib_name}
    PUBLIC $<BUILD_INTERFACE:${module_dir}> $<INSTALL_INTERFACE:include_${kind}>)

  # Handle OpenMP.
  if(OpenMP_Fortran_FOUND)
    target_link_libraries(${lib_name} PUBLIC OpenMP::OpenMP_Fortran)
  endif()

  # Link to sp library.
  target_link_libraries(${lib_name} PUBLIC sp::sp_${kind})
  
  list(APPEND LIB_TARGETS ${lib_name})
  
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
