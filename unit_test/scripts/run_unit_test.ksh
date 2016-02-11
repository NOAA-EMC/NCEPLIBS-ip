#!/bin/ksh

#---------------------------------------------------------------------------------
# Driver script to run the IPOLATES unit test.
#
# This script calls the "copygb" program to interpolate global datasets
# of scalar and vector data to several grids of various map projections
# using all ipolates interpolation options.  The output files (grib 1)
# are then compared to their corresponding "baseline" files using the
# "diffgb" utility.
#
# To run this script, type "./run.ksh"
#
# Three versions of copygb (located in ../src and ../exec) are tested:
#    - copygb_4 (uses single precision version of ipolates)
#    - copygb_8 (uses double precision version of ipolates)
#    - copygb_4 (uses mixed precision version of ipolates)
#
# The input data is located in the ../input_data directory.  There are two
# files (grib 1 format) so as to test the scalar and vector interpolation:
#    - ./scalar/global_snoalb.grb      - global snow albedo 
#    - ./vector/global_500mb_winds.grb - global u/v wind 
#
# The input data are interpolated to the following grids with the following
# ipolates interpolation options:
#    - NCEP grid 2 (global 2.5-deg lat/lon) using bilinear (IP option "0")
#    - NCEP grid 8 (mercator) using bicubic (IP option "1")
#    - NCEP grid 94 (rotated lat/lon "B") using neighbor (IP option "2")
#    - NCEP grid 98 (gaussian lat/lon) using budget (IP option "3")
#    - NCEP grid 100 (polar stereographic) using spectral (IP option "4")
#    - NCEP grid 192 (rotated lat/lon "E") using neighbor-budget (IP option "6")
#    - NCEP grid 211 (lambert conformal) using bilinear (IP option "0")
#
# The "baseline" data are located in subdirectories under ./baseline_data.
# The files are identified with the ncep grid number in the file name.
# The copygb program gives identical results for the double and mixed
# precision versions of ipolates.  These data are grib 1 format.  The 
# subdirectories are:
#    - ./scalar/4_byte: single precision albedo (scalar)
#    - ./scalar/8_byte: double and mixed precision albedo (scalar)
#    - ./vector/4_byte: single precision u/v wind (vector)
#    - ./vector/8_byte: double and mixed precision u/v wind (vector)
#
# The output from this script is piped to the screen.
#---------------------------------------------------------------------------------

#set -x

DIFFGB="../exec/diffgb.exe"

rm -f ../work/*.grb

for data_type in "scalar" "vector"  # scalar or vector test?
do
  case $data_type in  # the input data file
    "scalar")
      INPUT_FILE="../input_data/$data_type/global_snoalb.grb"  ;;
    "vector")
      INPUT_FILE="../input_data/$data_type/global_500mb_winds.grb"  ;;
  esac
  for precision in "4" "d" "8"  # test all three precision versions of ipolates
  do
    echo
    echo "***************************************************************"
    echo "*** TEST $precision BYTE VERSION OF IPOLATES LIBRARY FOR $data_type DATA ***"
    echo "***************************************************************"
    echo
    case $precision in  # the output files from copygb will be compared 
                        # to the baseline files in this directory.
      "4")
        BASELINE_DIR="../baseline_data/$data_type/4_byte"  ;;
      *)
        BASELINE_DIR="../baseline_data/$data_type/8_byte"  ;;
    esac
    COPYGB="../exec/copygb_$precision"
    for grid in 2 8 94 98 100 192 211  # the ncep grid number (defined above)
    do
      case $grid in   # the ipolates interpolation option (defined above)
        "2") 
          option="0" ;;
        "8")
          option="1" ;;
        "94")
          option="2" ;;
        "98")
          option="3" ;;
        "100")
          option="4 0 -1" ;;
        "192")
          option="6" ;;
        *)
          option="0" ;;
        esac
        echo
        echo "********************************************"
        echo "*** TEST NCEP GRID NUMBER $grid ***"
        echo "********************************************"
        echo
        $COPYGB -g$grid -i"${option}" -x $INPUT_FILE ../work/grid_${grid}.grb
        $DIFFGB -x ../work/grid_${grid}.grb $BASELINE_DIR/grid_${grid}.grb
        rm -f ../work/grid_${grid}.grb
    done
  done
done

exit 0
