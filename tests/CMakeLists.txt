# This is the CMake file for the test directory of NCEPLIBS-ip
#
# Kyle Gerheiser

#list(APPEND CMAKE_CTEST_ARGUMENTS "--output-on-failure")

# Copy test data to build directory.
file(COPY data/ DESTINATION data)

add_pfunit_ctest (interp_test
  TEST_SOURCES test_mod.pf
  LINK_LIBRARIES ip::ip_4 sp::sp_4
  )

# Run tests in the subdirectories.
#add_subdirectory(unit_test_for_c)
add_subdirectory(unit_test)
