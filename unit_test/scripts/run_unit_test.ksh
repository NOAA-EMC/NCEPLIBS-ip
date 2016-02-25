#!/bin/ksh

#---------------------------------------------------------------------------------
# Driver script to run the IPOLATES unit test.
#
# This script calls two programs to interpolate global datasets
# of scalar and vector data to several grids of various map projections
# using all ipolates interpolation options.  The interpolated data is
# then compared to "baseline" data, with a summary of differences 
# sent to standard output.
#
# To run this script interactively, type "./run_unit_test.ksh".
# On WCOSS, this script may be submitted to the compute nodes
# using "./run_wcoss.lsf".  On Theia, use "./run_theia.ksh"
#
# The program which tests the scalar interpolation is located in
# ../sorc/scalar.  After compilation, there are three executables
# in ../exec, one for each precision version of ipolates:
#    - scalar_4.exe (uses single precision version of ipolates)
#    - scalar_8.exe (uses double precision version of ipolates)
#    - scalar_d.exe (uses mixed precision version of ipolates)
#
# The program which tests the vector interpolation is located in
# ../sorc/vector.  As with the scalar program, after compilation
# there are three executables in ../exec:
#    - vector_4.exe (uses single precision version of ipolates)
#    - vector_8.exe (uses double precision version of ipolates)
#    - vector_d.exe (uses mixed precision version of ipolates)
#
# The input data is located in the ../input_data directory.  There are 
# two files (binary, little endian format), one with scalar data and one
# with vector data:
#    - ./scalar/global_snoalb.bin   - global snow albedo 
#    - ./vector/global_uv_wind.bin  - global 500mb u/v wind 
#
# The input data are interpolated to the following grids with the following
# ipolates interpolation options:
#    - grid 3 (global one-deg lat/lon) using bilinear (IP option "0")
#    - grid 8 (mercator) using bicubic (IP option "1")
#    - grid 127 (gaussian lat/lon) using neighbor (IP option "2")
#    - grid 203 (rotated lat/lon "E") using budget (IP option "3")
#    - grid 205 (rotated lat/lon "B") using spectral (IP option "4")
#    - grid 212 (polar stereographic) using neighbor-budget (IP option "6")
#    - grid 218 (lambert conformal) using bilinear (IP option "0")
#
# The "baseline" data are located in subdirectories under ./baseline_data.
# The files are identified with the grid number in the file name.
# Ipolates gives identical results for the double and mixed
# precision.  These data are binary little endian format.  The 
# subdirectories are:
#    - ./scalar/4_byte_bin: single precision albedo (scalar)
#    - ./scalar/8_byte_bin: double and mixed precision albedo (scalar)
#    - ./vector/4_byte_bin: single precision u/v wind (vector)
#    - ./vector/8_byte_bin: double and mixed precision u/v wind (vector)
#
# The output from this script is piped to the screen.
#---------------------------------------------------------------------------------

#set -x

cd ../work

ln -fs ../input_data/scalar/global_snoalb.bin   ./fort.9

for precision in "4" "d" "8"  # test all three precision versions of ipolates
do
  echo
  echo "***************************************************************"
  echo "*** TEST $precision BYTE VERSION OF IPOLATES LIBRARY FOR SCALAR DATA ***"
  echo "***************************************************************"

  EXEC="../exec/scalar_${precision}.exe"

  for grid in 3 8 127 203 205 212 218
  do
    case $grid in   # the ipolates interpolation option (defined above)
      "3")
        option="0" ;;
      "8")
        option="1" ;;
      "127")
        option="2" ;;
      "203")
        option="3" ;;
      "205")
        option="4" ;;
      "212")
        option="6" ;;
      *)
        option="0" ;;
    esac

    echo
    $EXEC $grid $option

  done
done

rm -f ./fort.9

ln -fs ../input_data/vector/global_uv_wind.bin   ./fort.9

for precision in "4" "d" "8"  # test all three precision versions of ipolates
do
  echo
  echo "***************************************************************"
  echo "*** TEST $precision BYTE VERSION OF IPOLATES LIBRARY FOR VECTOR DATA ***"
  echo "***************************************************************"

  EXEC="../exec/vector_${precision}.exe"

  for grid in 3 8 127 203 205 212 218
  do

    case $grid in   # the ipolates interpolation option (defined above)
      "3")
        option="0" ;;
      "8")
        option="1" ;;
      "127")
        option="2" ;;
      "203")
        option="3" ;;
      "205")
        option="4" ;;
      "212")
        option="6" ;;
      *)
        option="0" ;;
    esac

    echo
    $EXEC $grid $option

  done
done

rm -f ./fort.9

exit 0
