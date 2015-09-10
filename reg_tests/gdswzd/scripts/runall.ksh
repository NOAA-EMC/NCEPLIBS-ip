#!/bin/ksh

#--------------------------------------------------------------------
# Test the gdswzd suite of routines using a Fortran
# program.
#
# The program is compiled with all three byte versions
# of the 'control' and 'test' ip library.
#
# The three byte versions of the library are:
#  > 4 byte integer/4 byte float  ($bytesize=4)
#  > 8 byte integer/8 byte float  ($bytesize=8)
#  > 8 byte float/4 byte integer  ($bytesize=d)
#
# The gdswzd routines compute the following fields:
#
# - lat/lon from i/j OR i/j from lat lon
# - clockwise vector rotation sines/cosines
# - dx/dlon, dx/dlat, dy/dlon, dy/dlat 
# - grid box area 
#
# The routines are called twice for each grid to test both the
# i/j to lat/lon AND lat/lon to i/j transforms.  This is controled by setting
# the routine's IOPT argument to '0'/'-1'.  The transform should be reversable to
# within floating point differences.  If it is not reversable a warning
# message is printed to standard output and the regression test
# is considered failed.
#
# All fields computed from each call to gdswzd are output to a binary file.
# The file naming convention is: grid${gridnum}.iopt${0/m1}.bin.
# The files from the 'control' and 'test' ip libraries are compared
# and if not bit identical, the regression test fails.
#
# The grids tested are:
#
# 003        one-degree global lat/lon (ncep grid 3)
# 008        mercator (ncep grid 8)
# 127        t254 gaussian (ncep grid 127)
# 203h       rotated lat/lon e-staggered (number meaningless)
#            this is the old 12km eta grid - 'h' points
# 203v       rotated lat/lon e-staggered (number meaningless)
#            this is the old 12km eta grid - 'v' points
# 205h       rotated lat/lon b-staggered (number meaningless)
#            this is the 12km nam grid - 'h' points
# 205v       rotated lat/lon b-staggered (number meaningless)
#            this is the 12km nam grid - 'v' points
# 212        nh polar stereographic, spherical earth (number meaningless)
# 213        sh polar stereographic, spherical earth (number meaningless)
# 218        lambert conformal (ncep grid 218)
# 222        nh polar stereographic, elliptical earth (number meaningless)
#
# If any step in the regression test fails, the log file and
# any binary files will be stored in a "failed" sub-directory under
# working directory, $WORK_DIR.
#
# This script is run by the Runall.${machine}.ksh driver script located
# in ../reg_tests.
#--------------------------------------------------------------------

#set -x

echo
echo BEGIN GDSWZD REGRESSION TEST
echo

WORK_DIR=${WORK_DIR:-/stmpp1/$LOGNAME/regression}

REG_DIR=${REG_DIR:-../..}

# where the executables are located
EXEC_DIR=$REG_DIR/gdswzd/exec

WORK=${WORK_DIR}/gdswzd
rm -fr $WORK
mkdir -p $WORK
WORK_CTL=${WORK}/ctl
mkdir -p $WORK_CTL
cp $EXEC_DIR/gdswzd_ctl_*.exe  $WORK_CTL
WORK_TEST=${WORK}/test
mkdir -p $WORK_TEST
cp $EXEC_DIR/gdswzd_test_*.exe $WORK_TEST

reg_test_failed=0

for grids in "3" "8" "127" "203h" "203v" "212" "222" "213" "205h" "205v" "218"
do
  echo
  for bytesize in "4" "8" "d"  # test each library version
  do

    ctl_failed=0
    test_failed=0

    echo TEST GRID $grids ${bytesize}-BYTE FLOAT VERSION

    cd $WORK_CTL
    CTL_LOG=ctl.${bytesize}byte.grid${grids}.log
    gdswzd_ctl_${bytesize}.exe "$grids" > $CTL_LOG
    status=$?
# did 'control' executable run without error?
    if ((status != 0));then
      echo "** CONTROL RUN FAILED."
      ctl_failed=1
      reg_test_failed=1
    fi

    cd $WORK_TEST
    TEST_LOG=test.${bytesize}byte.grid${grids}.log
    gdswzd_test_${bytesize}.exe "$grids" > $TEST_LOG
    status=$?
# did 'test' executable run without error?
    if ((status != 0));then
      echo "** TEST RUN FAILED."
      test_failed=1
      reg_test_failed=1
    fi

# did the i/j to lat/lon transform fail for the control?
    grep -Eq 'BAD|ERROR' $WORK_CTL/$CTL_LOG
    status=$?
    if ((status == 0)); then
      echo "** PROBLEM WITH CTL RUN. CHECK LOG FILE."
      ctl_failed=1
      reg_test_failed=1
    fi

# did the i/j to lat/lon transform fail for the test?
    grep -Eq 'BAD|ERROR' $WORK_TEST/$TEST_LOG
    status=$?
    if ((status == 0)); then
      echo "** PROBLEM WITH TEST RUN. CHECK LOG FILE."
      test_failed=1
      reg_test_failed=1
    fi

# are binary files bit identical?
    if ((test_failed == 0 && ctl_failed == 0)); then
      cmp $WORK_CTL/grid${grids}.iopt0.bin $WORK_TEST/grid${grids}.iopt0.bin
      status=$?
      if ((status != 0));then
        echo "** BINARY FILES NOT BIT IDENTICAL. REGRESSION TEST FAILED."
        ctl_failed=1
        test_failed=1
        reg_test_failed=1
      fi
      cmp $WORK_CTL/grid${grids}.ioptm1.bin $WORK_TEST/grid${grids}.ioptm1.bin
      status=$?
      if ((status != 0));then
        echo "** BINARY FILES NOT BIT IDENTICAL. REGRESSION TEST FAILED."
        ctl_failed=1
        test_failed=1
        reg_test_failed=1
      fi
    fi

# if any step of regression test fails, save data in a subdirectory.

    if ((ctl_failed == 1));then
      FAILED_DIR=$WORK_CTL/failed.${bytesize}byte.grid${grids}
      mkdir -p $FAILED_DIR
      if [ -s $WORK_CTL/$CTL_LOG ]; then
        mv $WORK_CTL/$CTL_LOG $FAILED_DIR
      fi
      if [ -s $WORK_CTL/grid${grids}.iopt0.bin ];then
        mv $WORK_CTL/grid${grids}.iopt0.bin $FAILED_DIR
      fi
      if [ -s $WORK_CTL/grid${grids}.ioptm1.bin ];then
        mv $WORK_CTL/grid${grids}.ioptm1.bin $FAILED_DIR
      fi
      if [ -s $WORK_CTL/core* ];then
        mv $WORK_CTL/core*  $FAILED_DIR
      fi 
    fi

    if ((test_failed == 1));then
      FAILED_DIR=$WORK_TEST/failed.${bytesize}byte.grid${grids}
      mkdir -p $FAILED_DIR
      if [ -s $WORK_TEST/$TEST_LOG ];then
        mv $WORK_TEST/$TEST_LOG $FAILED_DIR
      fi
      if [ -s $WORK_TEST/grid${grids}.iopt0.bin ];then
        mv $WORK_TEST/grid${grids}.iopt0.bin $FAILED_DIR
      fi
      if [ -s $WORK_TEST/grid${grids}.ioptm1.bin ];then
        mv $WORK_TEST/grid${grids}.ioptm1.bin $FAILED_DIR
      fi
      if [ -s $WORK_TEST/core* ];then
        mv $WORK_TEST/core*  $FAILED_DIR
      fi
    fi

    rm -f $WORK_CTL/grid${grids}.iopt0.bin  $WORK_TEST/grid${grids}.iopt0.bin
    rm -f $WORK_CTL/grid${grids}.ioptm1.bin $WORK_TEST/grid${grids}.ioptm1.bin
    rm -f $WORK_CTL/$CTL_LOG $WORK_TEST/$TEST_LOG
    rm -f $WORK_CTL/core*  $WORK_TEST/core*

  done
done

if ((reg_test_failed == 0)); then
  echo
  echo "<<< GDSWZD REGRESSION TEST PASSED. >>>"
  echo
else
  echo
  echo "<<< GDSWZD REGRESSION TEST FAILED. >>>"
  echo
fi

exit 0
