# This is the CMake file for the test directory of NCEPLIBS-ip.
#
# Mark Potts, Kyle Gerheiser, Eric Engle

# Link data directory to find the test data.
execute_process(COMMAND cmake -E create_symlink
  "${CMAKE_CURRENT_SOURCE_DIR}/data"
  "${CMAKE_CURRENT_BINARY_DIR}/data" # New name
  )

# Set compiler flags for intel.
if(CMAKE_Fortran_COMPILER_ID MATCHES "^(Intel)$")
  set(CMAKE_Fortran_FLAGS "-r8 -g -traceback -check all -warn all -heap-arrays -assume byterecl ${CMAKE_Fortran_FLAGS} ")
  set(CMAKE_C_FLAGS "-std=c99")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "^(GNU)$")
  set(CMAKE_Fortran_FLAGS "-fdefault-real-8 -fno-range-check -g -fbacktrace -fcheck=all -Wall -O0 -fimplicit-none -Wsurprising -Wextra ${CMAKE_Fortran_FLAGS} ")
endif()

# Set compiler flags for GNU.
if(${CMAKE_Fortran_COMPILER_ID} MATCHES "^(GNU)$" AND ${CMAKE_Fortran_COMPILER_VERSION} VERSION_GREATER_EQUAL 10)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -w -fallow-argument-mismatch -fallow-invalid-boz")
endif()

foreach(kind "d")
  # Test earth_radius_mod.
  add_executable(test_earth_radius_${kind} test_earth_radius.F90)
  target_link_libraries(test_earth_radius_${kind} PUBLIC ip::ip_${kind})
  target_link_libraries(test_earth_radius_${kind} PUBLIC sp::sp_${kind})
  add_test(test_earth_radius_${kind} test_earth_radius_${kind})
  
  # grib-2 tests
  add_library(test_library_grib2_${kind} input_data_mod_grib2.F90 interp_mod_grib2.F90)
  target_link_libraries(test_library_grib2_${kind} PUBLIC ip::ip_${kind})
  target_link_libraries(test_library_grib2_${kind} PUBLIC sp::sp_${kind})
  
  add_executable(test_gdswzd_grib2_${kind} test_gdswzd_grib2.c)
  set_target_properties(test_gdswzd_grib2_${kind} PROPERTIES LINKER_LANGUAGE C)
  add_executable(test_scalar_grib2_${kind} test_scalar_grib2.F90)
  add_executable(test_vector_grib2_${kind} test_vector_grib2.F90)
  
  target_link_libraries(test_gdswzd_grib2_${kind} PRIVATE test_library_grib2_${kind})
  target_link_libraries(test_scalar_grib2_${kind} PRIVATE test_library_grib2_${kind})
  target_link_libraries(test_vector_grib2_${kind} PRIVATE test_library_grib2_${kind})
  
  add_test(test_gdswzd_c_grib2_${kind} test_gdswzd_grib2_${kind})
  
  # scalar tests
  add_test(test_lambert_bilinear_scalar_grib2_${kind} test_scalar_grib2_${kind} 218 0)
  add_test(test_gaussian_neighbor_scalar_grib2_${kind} test_scalar_grib2_${kind} 127 2)
  add_test(test_latlon_bilinear_scalar_grib2_${kind} test_scalar_grib2_${kind} 3 0)
  add_test(test_mercator_bicubic_scalar_grib2_${kind} test_scalar_grib2_${kind} 8 1)
  add_test(test_polar-stereo_neighbor-budget_scalar_grib2_${kind} test_scalar_grib2_${kind} 212 6)
  add_test(test_rotatedB_spectral_scalar_grib2_${kind} test_scalar_grib2_${kind} 205 4)
  add_test(test_rotatedE_budget_scalar_grib2_${kind} test_scalar_grib2_${kind} 203 3)

  # scalar station point tests
  add_test(test_station_points_bilinear_scalar_grib2_${kind} test_scalar_grib2_${kind} -1 0)
  add_test(test_station_points_bicubic_scalar_grib2_${kind} test_scalar_grib2_${kind} -1 1)
  add_test(test_station_points_neighbor_scalar_grib2_${kind} test_scalar_grib2_${kind} -1 2)
  add_test(test_station_points_budget_scalar_grib2_${kind} test_scalar_grib2_${kind} -1 3)
  add_test(test_station_points_spectral_scalar_grib2_${kind} test_scalar_grib2_${kind} -1 4)
#  add_test(test_station_points_neighbor-budget_scalar_grib2_${kind} test_scalar_grib2_${kind} -1 6)
  
  # # vector tests
  add_test(test_lambert_biliner_vector_grib2_${kind} test_vector_grib2_${kind} 218 0)
  add_test(test_gaussian_neighbor_vector_grib2_${kind} test_vector_grib2_${kind} 127 2)
  add_test(test_latlon_bilinear_vector_grib2_${kind} test_vector_grib2_${kind} 3 0)
  add_test(test_mercator_bicubic_vector_grib2_${kind} test_vector_grib2_${kind} 8 1)
  add_test(test_polar-stereo_neighbor-budget_vector_grib2_${kind} test_vector_grib2_${kind} 212 6)
  add_test(test_rotatedB_spectral_vector_grib2_${kind} test_vector_grib2_${kind} 205 4)
  add_test(test_rotatedE_budget_vector_grib2_${kind} test_vector_grib2_${kind} 203 3)

  # vector station point tests
  add_test(test_station_points_bilinear_vector_grib2_${kind} test_vector_grib2_${kind} -1 0)
  add_test(test_station_points_bicubic_vector_grib2_${kind} test_vector_grib2_${kind} -1 1)
  add_test(test_station_points_neighbor_vector_grib2_${kind} test_vector_grib2_${kind} -1 2)
  add_test(test_station_points_budget_vector_grib2_${kind} test_vector_grib2_${kind} -1 3)
  add_test(test_station_points_spectral_vector_grib2_${kind} test_vector_grib2_${kind} -1 4)
#  add_test(test_station_points_neighbor-budget_vector_grib2_${kind} test_vector_grib2_${kind} -1 6)
  
  # grib-1 tests
  add_library(test_library_grib1_${kind} input_data_mod_grib1.F90 interp_mod_grib1.F90)
  target_link_libraries(test_library_grib1_${kind} PUBLIC ip::ip_${kind})
  target_link_libraries(test_library_grib1_${kind} PUBLIC sp::sp_${kind})
  
  add_executable(test_gdswzd_grib1_${kind} test_gdswzd_grib1.c)
  set_target_properties(test_gdswzd_grib1_${kind} PROPERTIES LINKER_LANGUAGE C)
  add_executable(test_scalar_grib1_${kind} test_scalar_grib1.F90)
  add_executable(test_vector_grib1_${kind} test_vector_grib1.F90)
  
  target_link_libraries(test_scalar_grib1_${kind} PRIVATE test_library_grib1_${kind})
  target_link_libraries(test_vector_grib1_${kind} PRIVATE test_library_grib1_${kind})
  target_link_libraries(test_gdswzd_grib1_${kind} ip::ip_${kind})
  target_link_libraries(test_gdswzd_grib1_${kind} sp::sp_${kind})
  
  add_test(test_gdswzd_c_grib1_${kind} test_gdswzd_grib1_${kind})
  add_test(test_lambert_bilinear_scalar_grib1_${kind} test_scalar_grib1_${kind} 218 0)
  add_test(test_gaussian_neighbor_scalar_grib1_${kind} test_scalar_grib1_${kind} 127 2)
  add_test(test_latlon_bilinear_scalar_grib1_${kind} test_scalar_grib1_${kind} 3 0)
  add_test(test_mercator_bicubic_scalar_grib1_${kind} test_scalar_grib1_${kind} 8 1)
  add_test(test_polar-stereo_neighbor-budget_scalar_${kind} test_scalar_grib1_${kind} 212 6)
  add_test(test_rotatedB_spectral_scalar_${kind} test_scalar_grib1_${kind} 205 4)
  add_test(test_rotatedE_budget_scalar_grib1_${kind} test_scalar_grib1_${kind} 203 3)
  
  # vector tests
  add_test(test_lambert_biliner_vector_grib1_${kind} test_vector_grib1_${kind} 218 0)
  add_test(test_gaussian_neighbor_vector_grib1_${kind} test_vector_grib1_${kind} 127 2)
  add_test(test_latlon_bilinear_vector_grib1_${kind} test_vector_grib1_${kind} 3 0)
  add_test(test_mercator_bicubic_vector_grib1_${kind} test_vector_grib1_${kind} 8 1)
  add_test(test_polar-stereo_neighbor-budget_vector_grib1_${kind} test_vector_grib1_${kind} 212 6)
  add_test(test_rotatedB_spectral_vector_grib1_${kind} test_vector_grib1_${kind} 205 4)
  add_test(test_rotatedE_budget_vector_grib1_${kind} test_vector_grib1_${kind} 203 3)
endforeach()

