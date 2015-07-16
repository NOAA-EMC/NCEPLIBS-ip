#!/bin/ksh

#------------------------------------------------------------------------
# Run part 1 of the ipolates regression test to exercise the ipolates suite
# of routines.  These routines are exercised by a Fortran program.
#
# The program is compiled with all three byte versions
# of the 'control' and 'test' ip library.  The executables are located
# in the ./exec subdirectory (there are six).
#
# The three byte versions of the library are:
#  > 4 byte integer/4 byte float  ($bytesize=4)
#  > 8 byte integer/8 byte float  ($bytesize=8)
#  > 8 byte float/4 byte integer  ($bytesize=d)
#
# The program interpolates a global grid of substrate temperature to several 
# grids of various map projections.  The grids are:
#
#    3 - one-degree global lat/lon (ncep grid 3)
#    8 - mercator (ncep grid 8)
#  127 - t254 gaussian (ncep grid 127)
#  203 - rotated lat/lon e-staggered (number meaningless)
#        this is the old 12km eta grid - 'h' pts
#  205 - rotated lat/lon b-staggered (number meaningless)
#        this is the 12km nam grid - 'h' pts
#  212 - nh polar stereographic, spherical earth (number meaningless)
#  218 - lambert conformal (ncep grid 218)
#
# The input data is: ../data/global_tbot.1x1.grb2
# It is in grib 2 format.
#
# Use all possible ipolates interpolation options:
#
#    0 - bilinear
#    1 - bicubic
#    2 - neighbor
#    3 - budget
#    4 - spectral
#    6 - budget-neighbor
#
# This script is run by the Runall.${machine}.ksh driver
# script located in /reg_tests.
#
# The ipolates suite of routines contain threads.  Therefore, this 
# script is run twice, with 1 and 4 threads.  The number of 
# threads is passed in as an argument.  This is only used to
# name the work directory, and does not cause the regression test to
# run with that number of threads. The number of threads is set
# from the /reg_tests/Runall.${machine}.ksh driver script.  
#
# The interpolated temperature data is output to direct access
# binary files under WORK_DIR.  These files may be viewed in Grads
# using the control files in the ./grads subdirectory.  The file
# naming convention is:
#   grid${grid_num}.opt${interp_opt_num}.${bytesize}byte.bin
#
# Binary files of the interpolated data from the control and test
# iplibs are checked for bit-identicalness.  If not identical, the test
# is considered failed and the script will save the file in a
# working sub-directory with a "failed" extension.
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
  echo "BEGIN IPOLATES REGRESSION TEST WITH " $num_threads "THREADS"
else
  echo "ENTER NUMBER OF THREADS"
  exit 99
fi

WORK_DIR=${WORK_DIR:-/stmpp1/$LOGNAME/regression}

REG_DIR=${REG_DIR:-../..}

EXEC_DIR=$REG_DIR/ipolates/exec

INPUT_DATA=$REG_DIR/ipolates/data/global_tbot.1x1.grb2

WORK=$WORK_DIR/ipolates.${num_threads}threads
rm -fr $WORK
mkdir -p $WORK
WORK_CTL=${WORK}/ctl
mkdir -p $WORK_CTL
cp $EXEC_DIR/ipolates_ctl_*.exe  $WORK_CTL
cp $INPUT_DATA  $WORK_CTL/fort.9
WORK_TEST=${WORK}/test
mkdir -p $WORK_TEST
cp $EXEC_DIR/ipolates_test_*.exe $WORK_TEST
cp $INPUT_DATA  $WORK_TEST/fort.9

ulimit -s 2048000

reg_test_failed=0

for grids in "3" "8" "127" "203" "205" "212" "218" 
do
  echo
  for option in "0" "1" "2" "3" "4" "6"  # interpolation option
  do
    for bytesize in "4" "8" "d"   # bytesize version of library.
    do

      ctl_failed=0
      test_failed=0

      echo TEST ${bytesize}-BYTE VERSION FOR GRID $grids AND INTERP OPTION $option

      cd $WORK_CTL
      ipolates_ctl_${bytesize}.exe "$grids" "$option" > ctl.log
      status=$?
      if ((status != 0));then
        echo "** CONTROL RUN FAILED"
        ctl_failed=1
        reg_test_failed=1
      fi

      cd $WORK_TEST
      ipolates_test_${bytesize}.exe "$grids" "$option" > test.log
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
          echo "** BINARY FILES NOT BIT IDENTIAL."
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
  echo "<<< IPOLATES REGRESSION TEST WITH " $num_threads "THREADS PASSED. >>>"
  echo
else
  echo
  echo "<<< IPOLATES REGRESSION TEST WITH " $num_threads "THREADS FAILED. >>>"
  echo
fi

exit 0
