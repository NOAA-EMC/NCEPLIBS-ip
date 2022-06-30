# This is the CMake file for the test directory of NCEPLIBS-ip.
#
# Mark Potts, Kyle Gerheiser

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

# Test earth_radius_mod.
add_executable(test_earth_radius test_earth_radius.F90)
target_link_libraries(test_earth_radius PUBLIC ip::ip_d)
add_test(test_earth_radius test_earth_radius)

# grib-2 tests
add_library(test_library_grib2 input_data_mod_grib2.F90 interp_mod_grib2.F90)
target_link_libraries(test_library_grib2 PUBLIC ip::ip_d)

add_executable(test_gdswzd_grib2 test_gdswzd_grib2.c)
add_executable(test_scalar_grib2 test_scalar_grib2.F90)
add_executable(test_vector_grib2 test_vector_grib2.F90)

target_link_libraries(test_gdswzd_grib2 PRIVATE test_library_grib2)
target_link_libraries(test_scalar_grib2 PRIVATE test_library_grib2)
target_link_libraries(test_vector_grib2 PRIVATE test_library_grib2)

add_test(test_gdswzd_c_grib2 test_gdswzd_grib2)

# scalar tests
add_test(test_lambert_bilinear_scalar_grib2 test_scalar_grib2 218 0)
add_test(test_gaussian_neighbor_scalar_grib2 test_scalar_grib2 127 2)
add_test(test_latlon_bilinear_scalar_grib2 test_scalar_grib2 3 0)
add_test(test_mercator_bicubic_scalar_grib2 test_scalar_grib2 8 1)
add_test(test_polar-stereo_neighbor-budget_scalar_grib2 test_scalar_grib2 212 6)
add_test(test_rotatedB_spectral_scalar_grib2 test_scalar_grib2 205 4)
add_test(test_rotatedE_budget_scalar_grib2 test_scalar_grib2 203 3)

# # vector tests
add_test(test_lambert_biliner_vector_grib2 test_vector_grib2 218 0)
add_test(test_gaussian_neighbor_vector_grib2 test_vector_grib2 127 2)
add_test(test_latlon_bilinear_vector_grib2 test_vector_grib2 3 0)
add_test(test_mercator_bicubic_vector_grib2 test_vector_grib2 8 1)
add_test(test_polar-stereo_neighbor-budget_vector_grib2 test_vector_grib2 212 6)
add_test(test_rotatedB_spectral_vector_grib2 test_vector_grib2 205 4)
add_test(test_rotatedE_budget_vector_grib2 test_vector_grib2 203 3)


# grib-1 tests
add_library(test_library_grib1 input_data_mod_grib1.F90 interp_mod_grib1.F90)
target_link_libraries(test_library_grib1 PUBLIC ip::ip_d)

add_executable(test_gdswzd_grib1 test_gdswzd_grib1.c)
add_executable(test_scalar_grib1 test_scalar_grib1.F90)
add_executable(test_vector_grib1 test_vector_grib1.F90)

target_link_libraries(test_scalar_grib1 PRIVATE test_library_grib1)
target_link_libraries(test_vector_grib1 PRIVATE test_library_grib1)
target_link_libraries(test_gdswzd_grib1 ip::ip_d)

add_test(test_gdswzd_c_grib1 test_gdswzd_grib1)
add_test(test_lambert_bilinear_scalar_grib1 test_scalar_grib1 218 0)
add_test(test_gaussian_neighbor_scalar_grib1 test_scalar_grib1 127 2)
add_test(test_latlon_bilinear_scalar_grib1 test_scalar_grib1 3 0)
add_test(test_mercator_bicubic_scalar_grib1 test_scalar_grib1 8 1)
add_test(test_polar-stereo_neighbor-budget_scalar test_scalar_grib1 212 6)
add_test(test_rotatedB_spectral_scalar test_scalar_grib1 205 4)
add_test(test_rotatedE_budget_scalar_grib1 test_scalar_grib1 203 3)

# # vector tests
add_test(test_lambert_biliner_vector_grib1 test_vector_grib1 218 0)
add_test(test_gaussian_neighbor_vector_grib1 test_vector_grib1 127 2)
add_test(test_latlon_bilinear_vector_grib1 test_vector_grib1 3 0)
add_test(test_mercator_bicubic_vector_grib1 test_vector_grib1 8 1)
add_test(test_polar-stereo_neighbor-budget_vector_grib1 test_vector_grib1 212 6)
add_test(test_rotatedB_spectral_vector_grib1 test_vector_grib1 205 4)
add_test(test_rotatedE_budget_vector_grib1 test_vector_grib1 203 3)


