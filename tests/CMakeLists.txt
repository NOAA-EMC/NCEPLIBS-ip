# This is the CMake file for the test directory of NCEPLIBS-ip.
#
# Alex Richert, Mark Potts, Kyle Gerheiser, Eric Engle

function(create_sp_test name kind timeout)
  add_executable(${name}_${kind} ${name}.F90)

  # Include openMP if desired.
  if(OpenMP_Fortran_FOUND)
    target_link_libraries(${name}_${kind} PRIVATE OpenMP::OpenMP_Fortran)
  endif()
  target_link_libraries(${name}_${kind} PRIVATE ip::ip_${kind})
  if(CMAKE_Fortran_COMPILER_ID MATCHES "^(Intel|IntelLLVM)$")
    set_target_properties(${name}_${kind} PROPERTIES COMPILE_FLAGS "-convert big_endian ${fortran_${kind}_flags}")
  elseif(${CMAKE_Fortran_COMPILER_ID} MATCHES "^(GNU)$")
    set_target_properties(${name}_${kind} PROPERTIES COMPILE_FLAGS "-fconvert=big-endian ${fortran_${kind}_flags}")
  endif()
  add_test(NAME ${name}_${kind} COMMAND ${name}_${kind})
  target_compile_definitions(${name}_${kind} PRIVATE KIND_${kind})
  if(TEST_TIME_LIMIT)
    set_tests_properties(${name}_${kind} PROPERTIES TIMEOUT ${timeout})
  endif()
endfunction()

# Link data directory to find the test data.
execute_process(COMMAND cmake -E create_symlink
  "${CMAKE_CURRENT_SOURCE_DIR}/data"
  "${CMAKE_CURRENT_BINARY_DIR}/data" # New name
  )

# Set compiler flags.
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_DEBUG}")
if(CMAKE_Fortran_COMPILER_ID MATCHES "^(Intel|IntelLLVM)$")
  set(CMAKE_Fortran_FLAGS " -heap-arrays ${CMAKE_Fortran_FLAGS}")
  set(CMAKE_C_FLAGS "-std=c99")
endif()

foreach(kind ${kinds})
  string(TOUPPER ${kind} kind_definition)

  # Test ipxwafs routines
  add_executable(test_ipxwafs_${kind} test_ipxwafs.F90)
  target_link_libraries(test_ipxwafs_${kind} PUBLIC ip::ip_${kind})
  target_compile_definitions(test_ipxwafs_${kind} PRIVATE "LSIZE=${kind_definition}")
  set_target_properties(test_ipxwafs_${kind} PROPERTIES COMPILE_FLAGS "${fortran_${kind}_flags}")
  add_test(test_ipxwafs_${kind} test_ipxwafs_${kind})

  # Test earth_radius_mod.
  add_executable(test_earth_radius_${kind} test_earth_radius.F90)
  target_link_libraries(test_earth_radius_${kind} PUBLIC ip::ip_${kind})
  target_compile_definitions(test_earth_radius_${kind} PRIVATE "LSIZE=${kind_definition}")
  set_target_properties(test_earth_radius_${kind} PROPERTIES COMPILE_FLAGS "${fortran_${kind}_flags}")
  add_test(test_earth_radius_${kind} test_earth_radius_${kind})
  
  # grib-2 tests
  add_library(test_library_grib2_${kind} input_data_mod_grib2_${kind}.F90 interp_mod_grib2_${kind}.F90)
  target_link_libraries(test_library_grib2_${kind} PUBLIC ip::ip_${kind})
  target_compile_definitions(test_library_grib2_${kind} PRIVATE "LSIZE=${kind_definition}")
  set_target_properties(test_library_grib2_${kind} PROPERTIES COMPILE_FLAGS "${fortran_${kind}_flags}")
  
  add_executable(tst_gdswzd_grib2_${kind} tst_gdswzd_grib2.c)
  set_target_properties(tst_gdswzd_grib2_${kind} PROPERTIES LINKER_LANGUAGE C)
  target_compile_definitions(tst_gdswzd_grib2_${kind} PRIVATE "LSIZE=${kind_definition}")
  target_link_libraries(tst_gdswzd_grib2_${kind} PRIVATE test_library_grib2_${kind})

  add_executable(test_scalar_grib2_${kind} test_scalar_grib2.F90)
  add_executable(test_vector_grib2_${kind} test_vector_grib2.F90)
  target_link_libraries(test_scalar_grib2_${kind} PRIVATE test_library_grib2_${kind})
  target_link_libraries(test_vector_grib2_${kind} PRIVATE test_library_grib2_${kind})
  target_compile_definitions(test_scalar_grib2_${kind} PRIVATE "LSIZE=${kind_definition}")
  target_compile_definitions(test_vector_grib2_${kind} PRIVATE "LSIZE=${kind_definition}")
  set_target_properties(test_scalar_grib2_${kind} PROPERTIES COMPILE_FLAGS "${fortran_${kind}_flags}")
  set_target_properties(test_vector_grib2_${kind} PROPERTIES COMPILE_FLAGS "${fortran_${kind}_flags}")
  
  add_test(tst_gdswzd_c_grib2_${kind} tst_gdswzd_grib2_${kind})
  
  # scalar tests
  add_test(test_lambert_bilinear_scalar_grib2_${kind} test_scalar_grib2_${kind} 218 0)
  add_test(test_gaussian_neighbor_scalar_grib2_${kind} test_scalar_grib2_${kind} 127 2)
  add_test(test_latlon_bilinear_scalar_grib2_${kind} test_scalar_grib2_${kind} 3 0)
  add_test(test_mercator_bicubic_scalar_grib2_${kind} test_scalar_grib2_${kind} 8 1)
  add_test(test_polar_stereo_neighbor_budget_scalar_grib2_${kind} test_scalar_grib2_${kind} 212 6)
  add_test(test_rotatedB_spectral_scalar_grib2_${kind} test_scalar_grib2_${kind} 205 4)
  add_test(test_rotatedE_budget_scalar_grib2_${kind} test_scalar_grib2_${kind} 203 3)

  # scalar station point tests
  add_test(test_station_points_bilinear_scalar_grib2_${kind} test_scalar_grib2_${kind} -1 0)
  add_test(test_station_points_bicubic_scalar_grib2_${kind} test_scalar_grib2_${kind} -1 1)
  add_test(test_station_points_neighbor_scalar_grib2_${kind} test_scalar_grib2_${kind} -1 2)
  add_test(test_station_points_budget_scalar_grib2_${kind} test_scalar_grib2_${kind} -1 3)
  add_test(test_station_points_spectral_scalar_grib2_${kind} test_scalar_grib2_${kind} -1 4)
  add_test(test_station_points_neighbor_budget_scalar_grib2_${kind} test_scalar_grib2_${kind} -1 6)
  
  # # vector tests
  add_test(test_lambert_bilinear_vector_grib2_${kind} test_vector_grib2_${kind} 218 0)
  add_test(test_gaussian_neighbor_vector_grib2_${kind} test_vector_grib2_${kind} 127 2)
  add_test(test_latlon_bilinear_vector_grib2_${kind} test_vector_grib2_${kind} 3 0)
  add_test(test_mercator_bicubic_vector_grib2_${kind} test_vector_grib2_${kind} 8 1)
  add_test(test_polar_stereo_neighbor_budget_vector_grib2_${kind} test_vector_grib2_${kind} 212 6)
  add_test(test_rotatedB_spectral_vector_grib2_${kind} test_vector_grib2_${kind} 205 4)
  add_test(test_rotatedE_budget_vector_grib2_${kind} test_vector_grib2_${kind} 203 3)
  add_test(test_rotatedB_direct_spectral_vector_grib2_${kind} test_vector_grib2_${kind} 32769 4)
  add_test(test_rotatedB_direct_ncep_post_spectral_vector_grib2_${kind} test_vector_grib2_${kind} 32769b 4)
  add_test(test_rotatedE_direct_budget_vector_grib2_${kind} test_vector_grib2_${kind} 32768 3)

  # vector station point tests
  add_test(test_station_points_bilinear_vector_grib2_${kind} test_vector_grib2_${kind} -1 0)
  add_test(test_station_points_bicubic_vector_grib2_${kind} test_vector_grib2_${kind} -1 1)
  add_test(test_station_points_neighbor_vector_grib2_${kind} test_vector_grib2_${kind} -1 2)
  add_test(test_station_points_budget_vector_grib2_${kind} test_vector_grib2_${kind} -1 3)
  add_test(test_station_points_spectral_vector_grib2_${kind} test_vector_grib2_${kind} -1 4)
  add_test(test_station_points_neighbor_budget_vector_grib2_${kind} test_vector_grib2_${kind} -1 6)
  
  # grib-1 tests
  add_library(test_library_grib1_${kind} input_data_mod_grib1_${kind}.F90 interp_mod_grib1_${kind}.F90)
  target_link_libraries(test_library_grib1_${kind} PUBLIC ip::ip_${kind})
  target_compile_definitions(test_library_grib1_${kind} PRIVATE "LSIZE=${kind_definition}")
  set_target_properties(test_library_grib1_${kind} PROPERTIES COMPILE_FLAGS "${fortran_${kind}_flags}")
  
  add_executable(tst_gdswzd_grib1_${kind} tst_gdswzd_grib1.c)
  set_target_properties(tst_gdswzd_grib1_${kind} PROPERTIES LINKER_LANGUAGE C)
  target_compile_definitions(tst_gdswzd_grib1_${kind} PRIVATE "LSIZE=${kind_definition}")
  target_link_libraries(tst_gdswzd_grib1_${kind} ip::ip_${kind})

  add_executable(test_scalar_grib1_${kind} test_scalar_grib1.F90)
  add_executable(test_vector_grib1_${kind} test_vector_grib1.F90)
  target_link_libraries(test_scalar_grib1_${kind} PRIVATE test_library_grib1_${kind})
  target_link_libraries(test_vector_grib1_${kind} PRIVATE test_library_grib1_${kind})
  target_compile_definitions(test_scalar_grib1_${kind} PRIVATE "LSIZE=${kind_definition}")
  target_compile_definitions(test_vector_grib1_${kind} PRIVATE "LSIZE=${kind_definition}")
  set_target_properties(test_scalar_grib1_${kind} PROPERTIES COMPILE_FLAGS "${fortran_${kind}_flags}")
  set_target_properties(test_vector_grib1_${kind} PROPERTIES COMPILE_FLAGS "${fortran_${kind}_flags}")
  
  add_test(tst_gdswzd_c_grib1_${kind} tst_gdswzd_grib1_${kind})
  add_test(test_lambert_bilinear_scalar_grib1_${kind} test_scalar_grib1_${kind} 218 0)
  add_test(test_gaussian_neighbor_scalar_grib1_${kind} test_scalar_grib1_${kind} 127 2)
  add_test(test_latlon_bilinear_scalar_grib1_${kind} test_scalar_grib1_${kind} 3 0)
  add_test(test_mercator_bicubic_scalar_grib1_${kind} test_scalar_grib1_${kind} 8 1)
  add_test(test_polar_stereo_neighbor_budget_scalar_grib1_${kind} test_scalar_grib1_${kind} 212 6)
  add_test(test_rotatedB_spectral_scalar_grib1_${kind} test_scalar_grib1_${kind} 205 4)
  add_test(test_rotatedE_budget_scalar_grib1_${kind} test_scalar_grib1_${kind} 203 3)
  
  # vector tests
  add_test(test_lambert_bilinear_vector_grib1_${kind} test_vector_grib1_${kind} 218 0)
  add_test(test_gaussian_neighbor_vector_grib1_${kind} test_vector_grib1_${kind} 127 2)
  add_test(test_latlon_bilinear_vector_grib1_${kind} test_vector_grib1_${kind} 3 0)
  add_test(test_mercator_bicubic_vector_grib1_${kind} test_vector_grib1_${kind} 8 1)
  add_test(test_polar_stereo_neighbor_budget_vector_grib1_${kind} test_vector_grib1_${kind} 212 6)
  add_test(test_rotatedB_spectral_vector_grib1_${kind} test_vector_grib1_${kind} 205 4)
  add_test(test_rotatedE_budget_vector_grib1_${kind} test_vector_grib1_${kind} 203 3)

  # sp tests
  create_sp_test(test_ncpus ${kind} 0.3)
  create_sp_test(test_splaplac ${kind} 0.3)
  create_sp_test(test_splat ${kind} 0.3)
  create_sp_test(test_sppad ${kind} 0.3)
  create_sp_test(test_sptezv ${kind} 0.3)
  create_sp_test(test_fft ${kind} 0.3)
  create_sp_test(test_sptrung ${kind} 0.3)
  create_sp_test(test_sptrungv ${kind} 2)
endforeach()
