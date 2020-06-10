#list(APPEND CMAKE_CTEST_ARGUMENTS "--output-on-failure")

file(COPY data/ DESTINATION data)

find_package(sp REQUIRED)

add_pfunit_ctest (interp_test
  TEST_SOURCES test_mod.pf
  LINK_LIBRARIES ip_4 sp::sp_4
)
