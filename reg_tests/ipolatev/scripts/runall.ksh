#!/bin/ksh

#------------------------------------------------------------------------
# Run part 1 of the ipolatev regression test to exercise the ipolatev suite
# of routines.
#
# Interpolate a global grid of u/v wind to several
# grids of various map projections.  The grids are:
#
#    3 - one-degree global lat/lon (ncep grid 3)
#    8 - mercator (ncep grid 8)
#  127 - t254 gaussian (ncep grid 127)
#  203 - rotated lat/lon e-staggered (number refers to gds octet 6)
#  205 - rotated lat/lon b-staggered (number refers to gds octet 6)
#  212 - nh polar stereographic, spherical earth (number meaningless)
#  218 - lambert conformal (ncep grid 218)
#
#  Use all possible ipolatev interpolation options:
#
#    0 - bilinear
#    1 - bicubic
#    2 - neighbor
#    3 - budget
#    4 - spectral
#    6 - budget-neighbor
#
# The ipolatev suite of routines contain threads.  Therefore, this
# script is run twice, with 1 and 4 threads.  The number of
# threads is passed in as an argument.  This is only used to
# name the work directory, and does not cause the regression test to
# run with that number of threads. The number of threads is determined
# from the driver script.  You can run this script stand-alone.
# However, the default on ccs and zeus is to run with the maximum
# number of threads on a node.
#
# The control and test executables interpolate the data for
# a single grid and interpolation option.  Hence, they are invoked
# numerous times for all combinations of grids/interp options.
#
# The interpolated u/v wind data is output in a direct access
# binary file under WORK_DIR.  These files may be viewed in Grads
# using the control files in the ./grads subdirectory.  The file
# naming convention is "grid${grid_num}.opt${interp_opt_num}.bin"
#
# Binary files of the interpolated data from the control and test
# are check for bit-identicalness.  If not identical, the test
# is considered failed and the script will save the file in a
# working directory with a "failed" extension.
#
# This is part 1 of the regression test.  After completion of this
# script, the compare.ksh script is run to compare the binary
# files from the 1 and 4 thread runs.
#------------------------------------------------------------------------

#set -x

if (($# > 0))
then
  num_threads=$1
  echo
  echo "BEGIN IPOLATEV REGRESSION TEST WITH " $num_threads "THREADS"
else
  echo "ENTER NUMBER OF THREADS"
  exit 99
fi

WORK_DIR=${WORK_DIR:-/stmp/$LOGNAME/regression}

REG_DIR=${REG_DIR:-../..}

# where control and test programs are
EXEC_DIR=$REG_DIR/ipolatev/exec

# input u/v wind data
INPUT_DATA=$REG_DIR/ipolatev/data/gfs.500mb.winds.grb

WORK=${WORK_DIR}/ipolatev.${num_threads}threads
rm -fr $WORK
mkdir -p $WORK
WORK_CTL=${WORK}/ctl
mkdir -p $WORK_CTL
cp $EXEC_DIR/ctl/*.exe $WORK_CTL
cp $INPUT_DATA  $WORK_CTL/fort.9
WORK_TEST=${WORK}/test
mkdir -p $WORK_TEST
cp $EXEC_DIR/test/*.exe $WORK_TEST
cp $INPUT_DATA  $WORK_TEST/fort.9

ulimit -s 2048000

reg_test_failed=0

for grids in "3" "8" "127" "203" "205" "212" "218" 
do
  echo
  for option in "0" "1" "2" "3" "4" "6"  # interpolation option
  do
    for bytesize in "4" "8" "d"  # the three byte versions of the library
    do

      ctl_failed=0
      test_failed=0

      echo TEST ${bytesize}-BYTE VERSION FOR GRID $grids AND INTERP OPTION $option

      cd $WORK_CTL
      ipolatev_ctl_${bytesize}.exe "$grids" "$option" > ctl.log
      status=$?
      if ((status != 0));then
        echo "** CONTROL RUN FAILED"
        ctl_failed=1
        reg_test_failed=1
      fi

      cd $WORK_TEST
      ipolatev_test_${bytesize}.exe "$grids" "$option" > test.log
      status=$?
      if ((status != 0));then
        echo "** TEST RUN FAILED"
        test_failed=1
        reg_test_failed=1
      fi

      if ((test_failed == 0 && ctl_failed == 0)); then
        cmp $WORK_CTL/grid${grids}.opt${option}.bin $WORK_TEST/grid${grids}.opt${option}.bin
        status=$?
        if ((status != 0));then
          echo "** BINARY FILES NOT BIT IDENTIAL. TEST FAILED."
          ctl_failed=1
          test_failed=1
          reg_test_failed=1
        fi
      fi

      if ((ctl_failed == 1));then
        FAILED_DIR=$WORK_CTL/failed.grid${grids}.opt${option}.${bytesize}byte
        mkdir -p $FAILED_DIR
        if [ -s $WORK_CTL/ctl.log ]; then
          mv $WORK_CTL/ctl.log $FAILED_DIR
        fi
        if [ -s $WORK_CTL/grid${grids}.opt${option}.bin ];then
          mv $WORK_CTL/grid${grids}.opt${option}.bin  $FAILED_DIR
        fi
        if [ -s $WORK_CTL/core* ];then
          mv $WORK_CTL/core*  $FAILED_DIR
        fi
      else
        mv $WORK_CTL/grid${grids}.opt${option}.bin $WORK_CTL/grid${grids}.opt${option}.${bytesize}byte.bin
      fi

      if ((test_failed == 1));then
        FAILED_DIR=$WORK_TEST/failed.grid${grids}.opt${option}.${bytesize}byte
        mkdir -p $FAILED_DIR
        if [ -s $WORK_TEST/test.log ]; then
          mv $WORK_TEST/test.log $FAILED_DIR
        fi
        if [ -s $WORK_TEST/grid${grids}.opt${option}.bin ];then
          mv $WORK_TEST/grid${grids}.opt${option}.bin  $FAILED_DIR
        fi
        if [ -s $WORK_TEST/core* ];then
          mv $WORK_TEST/core*  $FAILED_DIR
        fi
      else
        mv $WORK_TEST/grid${grids}.opt${option}.bin $WORK_TEST/grid${grids}.opt${option}.${bytesize}byte.bin
      fi

      rm -f $WORK_CTL/ctl.log  $WORK_TEST/test.log

    done
  done
done

if ((reg_test_failed == 0));then
  echo
  echo "<<< IPOLATEV REGRESSION TEST WITH " $num_threads "THREADS PASSED. >>>"
  echo
else
  echo
  echo "<<< IPOLATEV REGRESSION TEST WITH " $num_threads "THREADS FAILED. >>>"
  echo
fi

exit 0
