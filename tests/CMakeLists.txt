# This is the CMake file for the tests/ directory of NCEPLIBS-ip.
# It define dozens of CTest unit tests. All new and modified code
# should have one or more unit tests to verify correct operation,
# including validation of outputs against known cases.
#
# Alex Richert, Mark Potts, Kyle Gerheiser, Eric Engle

if(FTP_TEST_FILES)
  set(TEST_TARBALL_FILE "ip-test-data-20241230.tgz")
  set(TEST_TARBALL_URL "https://ftp.emc.ncep.noaa.gov/static_files/public/NCEPLIBS-ip/${TEST_TARBALL_FILE}")
  if(EXISTS ${CMAKE_CURRENT_BINARY_DIR}/${TEST_TARBALL_FILE})
    message(STATUS "Using existing cached test file ${CMAKE_CURRENT_BINARY_DIR}/${TEST_TARBALL_FILE} ...")
  elseif(TEST_FILES_CACHE)
    if(EXISTS ${CMAKE_CURRENT_BINARY_DIR}/${TEST_TARBALL_FILE})
      message(STATUS "Using existing cached test file ${CMAKE_CURRENT_BINARY_DIR}/${TEST_TARBALL_FILE} ...")
    else()
      message(STATUS "Copying in cached test file ${TEST_FILES_CACHE} ...")
      file(COPY_FILE
        ${TEST_FILES_CACHE}
        ${CMAKE_CURRENT_BINARY_DIR}/${TEST_TARBALL_FILE}
      )
    endif()
  else()
    file(DOWNLOAD
      ${TEST_TARBALL_URL}
      ${CMAKE_CURRENT_BINARY_DIR}/${TEST_TARBALL_FILE}
      SHOW_PROGRESS
      STATUS DOWNLOAD_RESULT
      INACTIVITY_TIMEOUT 30
      )
    if(NOT DOWNLOAD_RESULT EQUAL 0)
      message(FATAL_ERROR "Failed to download ${TAR_FILE}")
    endif()
  endif()
  find_program(TAR tar)
  execute_process(
    COMMAND ${TAR} -xzf "${CMAKE_CURRENT_BINARY_DIR}/${TEST_TARBALL_FILE}" -C "${CMAKE_CURRENT_BINARY_DIR}"
    RESULT_VARIABLE TAR_RESULT
  )
  if(NOT TAR_RESULT EQUAL 0)
    message(FATAL_ERROR "Failed to extract ${TAR_FILE}")
  endif()
else()
  message(WARNING "FTP_TEST_FILES not enabled. Enable it or run 'ctest -L NO_INPUT_DATA' to avoid a bunch of test failures.")
endif()

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
  elseif(${CMAKE_Fortran_COMPILER_ID} MATCHES "^(LLVMFlang)$")
    set_target_properties(${name}_${kind} PROPERTIES COMPILE_FLAGS "-fconvert=big-endian ${fortran_${kind}_flags}")
  elseif(${CMAKE_Fortran_COMPILER_ID} MATCHES "^(NVHPC)$")
    set_target_properties(${name}_${kind} PROPERTIES COMPILE_FLAGS "-byteswapio ${fortran_${kind}_flags}")
  else()
    message(FATAL_ERROR "Unsupported Fortran compiler '${CMAKE_Fortran_COMPILER_ID}'")
  endif()
  add_test(NAME ${name}_${kind} COMMAND ${name}_${kind})
  target_compile_definitions(${name}_${kind} PRIVATE "LSIZE=${kind}")
  if(TEST_TIME_LIMIT)
    set_tests_properties(${name}_${kind} PROPERTIES TIMEOUT ${timeout})
  endif()
endfunction()

# Set compiler flags.
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_DEBUG}")
if(CMAKE_Fortran_COMPILER_ID MATCHES "^(Intel|IntelLLVM)$")
  set(CMAKE_Fortran_FLAGS " -heap-arrays ${CMAKE_Fortran_FLAGS}")
  set(CMAKE_C_FLAGS "-std=c99")
endif()

file(READ "${PROJECT_SOURCE_DIR}/VERSION" expected_version)
string(STRIP "${expected_version}" expected_version)

foreach(kind ${kinds})
  string(TOUPPER ${kind} kind_definition)

  ### Test version string interfaces
  add_executable(test_ip_version_c_${kind} test_ip_version.c)
  target_link_libraries(test_ip_version_c_${kind} PRIVATE ip::ip_${kind})
  add_test(NAME test_ip_version_c_${kind}
    COMMAND test_ip_version_c_${kind} "${expected_version}")

  add_executable(test_ip_version_fortran_${kind} test_ip_version.F90)
  target_link_libraries(test_ip_version_fortran_${kind} PRIVATE ip::ip_${kind})
  set_target_properties(test_ip_version_fortran_${kind} PROPERTIES COMPILE_FLAGS "${fortran_${kind}_flags}")
  add_test(NAME test_ip_version_fortran_${kind}
    COMMAND test_ip_version_fortran_${kind} "${expected_version}")

  ### Test ipxwafs routines
  add_executable(test_ipxwafs_${kind} test_ipxwafs.F90)
  target_link_libraries(test_ipxwafs_${kind} PUBLIC ip::ip_${kind})
  target_compile_definitions(test_ipxwafs_${kind} PRIVATE "LSIZE=${kind_definition}")
  set_target_properties(test_ipxwafs_${kind} PROPERTIES COMPILE_FLAGS "${fortran_${kind}_flags}")
  add_test(test_ipxwafs_${kind} test_ipxwafs_${kind})

  ### Test earth_radius_mod.
  add_executable(test_earth_radius_${kind} test_earth_radius.F90)
  target_link_libraries(test_earth_radius_${kind} PUBLIC ip::ip_${kind})
  target_compile_definitions(test_earth_radius_${kind} PRIVATE "LSIZE=${kind_definition}")
  set_target_properties(test_earth_radius_${kind} PROPERTIES COMPILE_FLAGS "${fortran_${kind}_flags}")
  add_test(test_earth_radius_${kind} test_earth_radius_${kind})

  ### GRIB2 tests
  # Create mod files by kind, GRIB1
  add_library(test_input_library_grib2_${kind} input_data_mod_grib2_${kind}.F90)
  target_link_libraries(test_input_library_grib2_${kind} PUBLIC ip::ip_${kind})
  target_compile_definitions(test_input_library_grib2_${kind} PUBLIC "LSIZE=${kind_definition}")
  set_target_properties(test_input_library_grib2_${kind} PROPERTIES COMPILE_FLAGS "${fortran_${kind}_flags}")
  add_library(test_library_grib2_${kind} interp_mod_grib2_${kind}.F90)
  target_link_libraries(test_library_grib2_${kind} PUBLIC test_input_library_grib2_${kind})
  set_target_properties(test_library_grib2_${kind} PROPERTIES COMPILE_FLAGS "${fortran_${kind}_flags}")

  # Set up executables, GRIB2
  add_executable(tst_gdswzd_grib2_${kind} tst_gdswzd_grib2.c)
  set_target_properties(tst_gdswzd_grib2_${kind} PROPERTIES LINKER_LANGUAGE C)
  target_compile_definitions(tst_gdswzd_grib2_${kind} PRIVATE "LSIZE=${kind_definition}")
  target_link_libraries(tst_gdswzd_grib2_${kind} PRIVATE test_library_grib2_${kind})

  add_executable(test_scalar_grib2_${kind} test_scalar_grib2.F90)
  add_executable(test_vector_grib2_${kind} test_vector_grib2.F90)
  target_link_libraries(test_scalar_grib2_${kind} PRIVATE test_library_grib2_${kind})
  target_link_libraries(test_vector_grib2_${kind} PRIVATE test_library_grib2_${kind})
  set_target_properties(test_scalar_grib2_${kind} PROPERTIES COMPILE_FLAGS "${fortran_${kind}_flags}")
  set_target_properties(test_vector_grib2_${kind} PROPERTIES COMPILE_FLAGS "${fortran_${kind}_flags}")

  # Test for ncep_post_arakawa C interface
  add_executable(tst_ncep_post_arakawa_${kind} tst_ncep_post_arakawa.c)
  target_link_libraries(tst_ncep_post_arakawa_${kind} PUBLIC ip::ip_${kind})
  set_target_properties(tst_ncep_post_arakawa_${kind} PROPERTIES LINKER_LANGUAGE C)
  target_compile_definitions(tst_ncep_post_arakawa_${kind} PRIVATE "LSIZE=${kind_definition}")
  add_test(tst_ncep_post_arakawa_${kind} tst_ncep_post_arakawa_${kind})

  # GDSWZD test, GRIB2
  add_test(tst_gdswzd_c_grib2_${kind} tst_gdswzd_grib2_${kind})

  # Scalar tests, GRIB2
  add_test(test_lambert_bilinear_scalar_grib2_${kind} test_scalar_grib2_${kind} 218 0)
  add_test(test_gaussian_neighbor_scalar_grib2_${kind} test_scalar_grib2_${kind} 127 2)
  add_test(test_latlon_bilinear_scalar_grib2_${kind} test_scalar_grib2_${kind} 3 0)
  add_test(test_mercator_bicubic_scalar_grib2_${kind} test_scalar_grib2_${kind} 8 1)
  add_test(test_polar_stereo_neighbor_budget_scalar_grib2_${kind} test_scalar_grib2_${kind} 212 6)
  add_test(test_rotatedB_spectral_scalar_grib2_${kind} test_scalar_grib2_${kind} 205 4)
  add_test(test_rotatedE_budget_scalar_grib2_${kind} test_scalar_grib2_${kind} 203 3)

  # Scalar station point tests, GRIB2
  add_test(test_station_points_bilinear_scalar_grib2_${kind} test_scalar_grib2_${kind} -1 0)
  add_test(test_station_points_bicubic_scalar_grib2_${kind} test_scalar_grib2_${kind} -1 1)
  add_test(test_station_points_neighbor_scalar_grib2_${kind} test_scalar_grib2_${kind} -1 2)
  add_test(test_station_points_budget_scalar_grib2_${kind} test_scalar_grib2_${kind} -1 3)
  add_test(test_station_points_spectral_scalar_grib2_${kind} test_scalar_grib2_${kind} -1 4)
  add_test(test_station_points_neighbor_budget_scalar_grib2_${kind} test_scalar_grib2_${kind} -1 6)

  # Vector tests, GRIB2
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

  # Vector station point tests, GRIB2
  add_test(test_station_points_bilinear_vector_grib2_${kind} test_vector_grib2_${kind} -1 0)
  add_test(test_station_points_bicubic_vector_grib2_${kind} test_vector_grib2_${kind} -1 1)
  add_test(test_station_points_neighbor_vector_grib2_${kind} test_vector_grib2_${kind} -1 2)
  add_test(test_station_points_budget_vector_grib2_${kind} test_vector_grib2_${kind} -1 3)
  add_test(test_station_points_spectral_vector_grib2_${kind} test_vector_grib2_${kind} -1 4)
  add_test(test_station_points_neighbor_budget_vector_grib2_${kind} test_vector_grib2_${kind} -1 6)

  ### GRIB1 tests
  # Create mod files by kind, GRIB1
  add_library(test_input_library_grib1_${kind} input_data_mod_grib1_${kind}.F90)
  target_link_libraries(test_input_library_grib1_${kind} PUBLIC ip::ip_${kind})
  target_compile_definitions(test_input_library_grib1_${kind} PUBLIC "LSIZE=${kind_definition}")
  set_target_properties(test_input_library_grib1_${kind} PROPERTIES COMPILE_FLAGS "${fortran_${kind}_flags}")
  add_library(test_library_grib1_${kind} interp_mod_grib1_${kind}.F90)
  target_link_libraries(test_library_grib1_${kind} PUBLIC test_input_library_grib1_${kind})
  set_target_properties(test_library_grib1_${kind} PROPERTIES COMPILE_FLAGS "${fortran_${kind}_flags}")

  # Set up executables, GRIB1
  add_executable(tst_gdswzd_grib1_${kind} tst_gdswzd_grib1.c)
  set_target_properties(tst_gdswzd_grib1_${kind} PROPERTIES LINKER_LANGUAGE C)
  target_compile_definitions(tst_gdswzd_grib1_${kind} PRIVATE "LSIZE=${kind_definition}")
  target_link_libraries(tst_gdswzd_grib1_${kind} ip::ip_${kind})

  add_executable(test_scalar_grib1_${kind} test_scalar_grib1.F90)
  add_executable(test_vector_grib1_${kind} test_vector_grib1.F90)
  target_link_libraries(test_scalar_grib1_${kind} PRIVATE test_library_grib1_${kind})
  target_link_libraries(test_vector_grib1_${kind} PRIVATE test_library_grib1_${kind})
  set_target_properties(test_scalar_grib1_${kind} PROPERTIES COMPILE_FLAGS "${fortran_${kind}_flags}")
  set_target_properties(test_vector_grib1_${kind} PROPERTIES COMPILE_FLAGS "${fortran_${kind}_flags}")

  # GDSWZD test, GRIB1
  add_test(tst_gdswzd_c_grib1_${kind} tst_gdswzd_grib1_${kind})

  # Scalar tests, GRIB1
  add_test(test_lambert_bilinear_scalar_grib1_${kind} test_scalar_grib1_${kind} 218 0)
  add_test(test_gaussian_neighbor_scalar_grib1_${kind} test_scalar_grib1_${kind} 127 2)
  add_test(test_latlon_bilinear_scalar_grib1_${kind} test_scalar_grib1_${kind} 3 0)
  add_test(test_mercator_bicubic_scalar_grib1_${kind} test_scalar_grib1_${kind} 8 1)
  add_test(test_polar_stereo_neighbor_budget_scalar_grib1_${kind} test_scalar_grib1_${kind} 212 6)
  add_test(test_rotatedB_spectral_scalar_grib1_${kind} test_scalar_grib1_${kind} 205 4)
  add_test(test_rotatedE_budget_scalar_grib1_${kind} test_scalar_grib1_${kind} 203 3)

  # Vector tests, GRIB1
  add_test(test_lambert_bilinear_vector_grib1_${kind} test_vector_grib1_${kind} 218 0)
  add_test(test_gaussian_neighbor_vector_grib1_${kind} test_vector_grib1_${kind} 127 2)
  add_test(test_latlon_bilinear_vector_grib1_${kind} test_vector_grib1_${kind} 3 0)
  add_test(test_mercator_bicubic_vector_grib1_${kind} test_vector_grib1_${kind} 8 1)
  add_test(test_polar_stereo_neighbor_budget_vector_grib1_${kind} test_vector_grib1_${kind} 212 6)
  add_test(test_rotatedB_spectral_vector_grib1_${kind} test_vector_grib1_${kind} 205 4)
  add_test(test_rotatedE_budget_vector_grib1_${kind} test_vector_grib1_${kind} 203 3)

  ### Various sp_mod tests
  create_sp_test(test_ncpus ${kind} 0.3)
  create_sp_test(test_splaplac ${kind} 0.3)
  create_sp_test(test_splat ${kind} 0.3)
  create_sp_test(test_sppad ${kind} 0.3)
  create_sp_test(test_sptezv ${kind} 0.3)
  create_sp_test(test_fft ${kind} 0.3)
  create_sp_test(test_sptrung ${kind} 0.3)
  create_sp_test(test_sptrungv ${kind} 2)

  ### Tests for fftpack subroutines
  add_library(test_fftpack_mod_${kind} test_fftpack_mod_${kind}.F90)
  target_link_libraries(test_fftpack_mod_${kind} PRIVATE ip::ip_${kind})
  set_target_properties(test_fftpack_mod_${kind} PROPERTIES COMPILE_FLAGS "${fortran_${kind}_flags}")
  add_executable(test_fftpack_${kind} test_fftpack.F90)
  target_compile_definitions(test_fftpack_${kind} PRIVATE "LSIZE=${kind}")
  set_target_properties(test_fftpack_${kind} PROPERTIES COMPILE_FLAGS "${fortran_${kind}_flags}")
  target_link_libraries(test_fftpack_${kind} PRIVATE test_fftpack_mod_${kind})
  target_include_directories(test_fftpack_${kind} PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
  target_include_directories(test_fftpack_${kind} PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/../src/include_${kind})
  add_test(test_fftpack_${kind} test_fftpack_${kind})

  ### Add test labels:

  # Label slow-running tests so they can be skipped using 'ctest -LE/--label-exclude'
  set(slowtests
    "test_rotatedB_spectral_vector_grib1"
    "test_rotatedB_spectral_scalar_grib1"
    "test_rotatedB_spectral_vector_grib2"
    "test_rotatedB_direct_spectral_vector_grib2"
    "test_rotatedB_direct_ncep_post_spectral_vector_grib2"
    "test_rotatedB_spectral_scalar_grib2"
  )

  foreach(slowtest IN LISTS slowtests)
    set_tests_properties("${slowtest}_${kind}" PROPERTIES LABELS "SLOW_TEST")
  endforeach()

  # Label tests without input data to allow trimmed down unit test runs
  # using 'ctest -L NO_INPUT_DATA'.
  set(noinputdatatests
    "test_earth_radius"
    "test_fft"
    "test_ip_version_c"
    "test_ip_version_fortran"
    "test_ipxwafs"
    "test_ncpus"
    "test_splaplac"
    "test_splat"
    "test_sppad"
    "test_sptezv"
    "tst_gdswzd_c_grib1"
    "tst_gdswzd_c_grib2"
    "tst_ncep_post_arakawa"
  )
  foreach(noinputdatatest IN LISTS noinputdatatests)
    set_tests_properties("${noinputdatatest}_${kind}" PROPERTIES LABELS "NO_INPUT_DATA")
  endforeach()

endforeach()
