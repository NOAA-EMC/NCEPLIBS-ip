#!/bin/ksh

#--------------------------------------------------------------------
# Test the gdswiz and gdswzd suite of routines.
#
# Binary files from the 'control' and 'test' containing gdswiz
# and gdswwzd related fields, such as lat/lon, are checked for
# bit-identicalness.  If not identical, the regression test is 
# considered failed.  Each run outputs two binary files -
# one with fields computed with the gdswiz/wzd iopt option
# set to '0', and one with the iopt option set to '-1'.
# (The files contain "iopt0/m1" in their name.)
# See comments in the source code for more details.
#
# The i/j to lat/lon transform should be reversable.
# If not, an error message is printed to the text log file.
# When this happens, the regression test is considered failed.
#
# If any step in the regression test fails, the log file and
# any binary files will be stored in a "failed" sub-directory under
# $WORK_DIR.
#--------------------------------------------------------------------

#set -x

echo
echo BEGIN GDSWIZ/WZD REGRESSION TEST
echo

WORK_DIR=${WORK_DIR:-/stmp/$LOGNAME/regression}

REG_DIR=${REG_DIR:-../..}

# where the executables are located
EXEC_DIR=$REG_DIR/gdswiz_wzd/exec

WORK=${WORK_DIR}/gdswiz_wzd
rm -fr $WORK
mkdir -p $WORK
WORK_CTL=${WORK}/ctl
mkdir -p $WORK_CTL
cp $EXEC_DIR/gdswiz_wzd_ctl_*.exe  $WORK_CTL
WORK_TEST=${WORK}/test
mkdir -p $WORK_TEST
cp $EXEC_DIR/gdswiz_wzd_test_*.exe $WORK_TEST

reg_test_failed=0

for routine in "WIZ" "WZD"  # test gdswiz and gdswzd separately
do
  echo
  echo RUN REGRESSION TEST FOR GDS${routine} ROUTINES 
  for grids in "3" "8" "203" "127" "212" "213" "218" "205" "201" "202" "222"
  do
    echo
    for bytesize in "4" "8" "d"  # test each library version
    do

      ctl_failed=0
      test_failed=0

      echo TEST GRID $grids ${bytesize}-BYTE FLOAT VERSION

      cd $WORK_CTL
      CTL_LOG=ctl.${routine}.${bytesize}byte.grid${grids}.log
      gdswiz_wzd_ctl_${bytesize}.exe "$routine" "$grids" > $CTL_LOG
      status=$?
# did 'control' executable run without error?
      if ((status != 0));then
        echo "** CONTROL RUN FAILED."
        ctl_failed=1
        reg_test_failed=1
      fi

      cd $WORK_TEST
      TEST_LOG=test.${routine}.${bytesize}byte.grid${grids}.log
      gdswiz_wzd_test_${bytesize}.exe "$routine" "$grids" > $TEST_LOG
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
        FAILED_DIR=$WORK_CTL/failed.${routine}.${bytesize}byte.grid${grids}
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
        FAILED_DIR=$WORK_TEST/failed.${routine}.${bytesize}byte.grid${grids}
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
done

if ((reg_test_failed == 0)); then
  echo
  echo "<<< GDSWIZ/WZD REGRESSION TEST PASSED. >>>"
  echo
else
  echo
  echo "<<< GDSWIZ/WZD REGRESSION TEST FAILED. >>>"
  echo
fi

exit 0
