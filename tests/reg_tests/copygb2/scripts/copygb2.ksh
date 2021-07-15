#!/bin/ksh

#--------------------------------------------------------------------
# Invoke a modified version of copygb2 to interpolate global
# lat/lon scalar and vector data to several grids of different
# map projections.  copygb2 invokes routine "ipolates" for scalar
# interpolation and "ipolatev" for vector interpolation.  All
# interpolation options are used:
#    - bilinear
#    - bicubic 
#    - neighbor
#    - budget  
#    - spectral 
#    - neighbor-budget 
#
# The output grids are:
#    - regular lat/lon
#    - mercator
#    - t254 gaussian 
#    - nh polar stereographic 
#    - lambert conformal
#
# The input data - in grib 2 - are located in the ../data directory:
#    - elevation.grb2  (global 1-deg terrain height)
#    - mxsnoalb.grb2   (global 1-deg max snow albedo - land only)
#    - uv_wind.grb2    (global 1-deg 500 mb u/v wind)
#
# The vector interpolation uses the u/v wind file as input.
# The scalar interpolation uses the terrain height (spectral
# interpolation option) or maximum snow albedo (all other options)
# as input.  The spectral interpolation does not work with
# land only (i.e., bitmapped) data.
#
# copygb2 is compiled with all three byte versions
# of the 'control' and 'test' ip2 library:
#  > 4 byte integer/4 byte float  
#  > 8 byte integer/8 byte float 
#  > 8 byte float/4 byte integer  
#
# The executables (a total of six) are stored in the
# ./exec subdirectory.
#
# The interpolated output from copygb2 is a grib 2 file.
# If the output grib 2 files from the 'test' and 'control' ip2libs
# are not bit identical, the regression test fails.
#
# Some script variable defintions:
#   $bytesize - ip2lib byte precision: 
#               "4" is 4 byte integer/4 byte float
#               "8" is 8 byte integer/8 byte float
#               "d" is 8 byte float/4 byte integer
#   $gridnum  - "4"   regular 0.5-deg lat/lon
#               "8"   mercator
#               "127" t254 gaussian 
#               "212" nh polar stereographic 
#               "218" lambert conformal
#   $int_type - "scalar"; perform scalar interpolation
#               "vector"; perform vector interpolation
#   $option   - interpolation option:
#               "0" bilinear
#               "1" bicubic
#               "2" neighbor
#               "3" budget
#               "4" spectral
#               "6" neighbor-budget
#               
# This is considered a 'supplemental' regression test and is not
# run as part of the full suite of tests from the main
# /reg_tests/Runall.${machine}.ksh driver scripts.  Instead,
# it is run stand-alone using the driver scripts in the
# /copygb2/scripts directory.
#--------------------------------------------------------------------

#set -x

echo
echo "BEGIN COPYGB2 REGRESSION TEST"

APRUN=${APRUN:-" "}

REG_DIR=${REG_DIR:-$(pwd)/../..}
EXEC_DIR=$REG_DIR/copygb2/exec
DATA_DIR=$REG_DIR/copygb2/data

WORK_DIR=${WORK_DIR:-/stmpp1/$LOGNAME/regression}
WORK=${WORK_DIR}/copygb2
rm -fr $WORK
mkdir -p $WORK
WORK_CTL=${WORK}/ctl
mkdir -p $WORK_CTL
WORK_TEST=${WORK}/test
mkdir -p $WORK_TEST

cp $EXEC_DIR/copygb2_ctl_*  $WORK_CTL
cp $EXEC_DIR/copygb2_test_* $WORK_TEST

# The output grids.  These are defined by their grid definition templates.
# ncep grid 4 - regular lat/lon
grid[4]="0 6 255 -1 255 -1 255 -1 720 361 0 -1 90000000 0 56 -90000000 359500000 500000 500000 0"

# ncep grid 8 - mercator
grid[8]="10 6 255 -1 255 -1 255 -1 116 44 -48670000 3104000 56 22500000 61050000 0 64 0 318830000 318830000"

# ncep grid 127 - t254 gaussain
grid[127]="40 6 255 -1 255 -1 255 -1 768 384 0 -1 89642000 0 48 -89642000 359531000 469000 192 0"

# afwa grid 212 - nh polar stereographic
grid[212]="20 6 255 -1 255 -1 255 -1 512 512 -20826000 145000000 56 60000000 280000000 47625000 47625000 0 0"

# ncep grid 218 - lambert conformal
grid[218]="30 6 255 -1 255 -1 255 -1 614 428 12190000 226541000 56 25000000 265000000 12191000 12191000 0 64 25000000 25000000 -90000000 0"

reg_test_failed=0

for bytesize in 4 d 8          # test all three byte versions of library.
do
  echo
  echo TEST $bytesize BYTE VERSION OF LIBRARY.
  for int_type in scalar vector   # test scalar and vector data
  do
    echo
    echo INTERPOLATE $int_type DATA.
    echo
    for option in 0 1 2 3 "4 0 -1" 6  # interpolation option
    do
      if [[ $option == "4 0 -1" ]]; then # spectral option does not work with bitmaps.
        if [[ $int_type == "scalar" ]]; then
          INPUT_FILE=$DATA_DIR/elevation.grb2
        else
          INPUT_FILE=$DATA_DIR/uv_wind.grb2
        fi
        opt=4
      else    # non-spectral option
        if [[ $int_type == "scalar" ]]; then
          INPUT_FILE=$DATA_DIR/mxsnoalb.grb2
        else
          INPUT_FILE=$DATA_DIR/uv_wind.grb2
        fi 
        opt=$option
      fi
      for gridnum in 4 8 127 212 218  # output grid
      do
        echo PROCESS GRID ${gridnum} USING INTERPOLATION OPTION ${opt}
        ctl_failed=0
        test_failed=0
        OUTPUT_FILE=grid${gridnum}.opt${opt}.${int_type}.${bytesize}byte.grb2
        cd $WORK_CTL
        $APRUN ./copygb2_ctl_${bytesize} -g "${grid[gridnum]}" -i"${option}" -x $INPUT_FILE $OUTPUT_FILE
        status=$?
        if ((status != 0)); then
          echo "** PROBLEM WITH CONTROL RUN."
          reg_test_failed=1
          ctl_failed=1
        fi
        if [[ ! -s $OUTPUT_FILE ]];then
          echo "** CONTROL RUN CREATED ZERO BYTE FILE"
          reg_test_failed=1
          ctl_failed=1
        fi
        cd $WORK_TEST
        $APRUN ./copygb2_test_${bytesize} -g "${grid[gridnum]}" -i"${option}" -x $INPUT_FILE $OUTPUT_FILE
        status=$?
        if ((status != 0)); then
          echo "** PROBLEM WITH TEST RUN."
          reg_test_failed=1
          test_failed=1
        fi
        if [[ ! -s $OUTPUT_FILE ]];then
          echo "** TEST RUN CREATED ZERO BYTE FILE"
          reg_test_failed=1
          test_failed=1
        fi
        cd $WORK
        if ((ctl_failed == 0 && test_failed == 0));then
          cmp $WORK_CTL/$OUTPUT_FILE  $WORK_TEST/$OUTPUT_FILE
          status=$?
          if ((status != 0)); then
            echo "** GRIB FILES NOT BIT IDENTICAL. REGRESSION TEST FAILS."
            mv $WORK_CTL/$OUTPUT_FILE $WORK_CTL/${OUTPUT_FILE}.failed
            mv $WORK_TEST/$OUTPUT_FILE $WORK_TEST/${OUTPUT_FILE}.failed
            reg_test_failed=1
          fi
        fi
        rm -f $WORK_CTL/$OUTPUT_FILE
        rm -f $WORK_TEST/$OUTPUT_FILE
      done
    done
  done
done

if ((reg_test_failed == 0)); then
  echo
  echo "<<< COPYGB2 REGRESSION TEST PASSED. >>>"
  echo
else
  echo
  echo "<<< COPYGB2 REGRESSION TEST FAILED. >>>"
  echo
fi

exit 0
