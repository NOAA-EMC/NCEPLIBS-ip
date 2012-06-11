#!/bin/ksh

#--------------------------------------------------------------------
# Test the gdswiz and gdswzd suite of routines.
#
# Binary files from the 'control' and 'test' containing gdswiz
# and gdswwzd related fields, such as lat/lon, are checked for
# bit-identicalness.  If not identical, the regression test is 
# considered failed.
#
# The i/j to lat/lon transform should be reversable.
# If not, an error message is printed to the text log file.
# When this happens, the regression test is considered failed.
#
# If any step in the regression test fails, the log file and
# any binary files will be stored in a "failed" directory under
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
cp $EXEC_DIR/ctl/*exe $WORK_CTL
WORK_TEST=${WORK}/test
mkdir -p $WORK_TEST
cp $EXEC_DIR/test/*exe $WORK_TEST

failed=0

for routine in "WIZ" "WZD"  # test gdswiz and gdswzd separately
do
  echo
  echo RUN REGRESSION TEST FOR GDS${routine} ROUTINES 
  for grids in "3" "8" "203" "127" "212" "213" "218" "205" "201" "202" "222"
  do
    echo
    for bytesize in "4" "8" "d"  # test each library version
    do

      save_ctl=0
      save_test=0

      echo TEST GRID $grids ${bytesize}-BYTE FLOAT VERSION

      cd $WORK_CTL
      CTL_LOG=ctl.${routine}.${bytesize}byte.grid${grids}.log
      gdswiz_wzd_ctl_${bytesize}.exe "$routine" "$grids" > $CTL_LOG
# did 'control' executable run without error?
      status=$?
      if ((status != 0));then
        echo "** CONTROL RUN FAILED."
        save_ctl=1
        failed=1
      fi

      cd $WORK_TEST
      TEST_LOG=test.${routine}.${bytesize}byte.grid${grids}.log
      gdswiz_wzd_test_${bytesize}.exe "$routine" "$grids" > $TEST_LOG
# did 'test' executable run without error?
      if ((status != 0));then
        echo "** TEST RUN FAILED."
        save_test=1
        failed=1
      fi

# are binary files bit identical?
      cmp $WORK_CTL/grid${grids}.bin $WORK_TEST/grid${grids}.bin
      status=$?
      if ((status != 0));then
        echo "** BINARY FILES NOT BIT IDENTICAL. TEST FAILED."
        save_ctl=1
        save_test=1
        failed=1
      fi

# did the i/j to lat/lon transform fail for the control?
      grep -Eq 'BAD|ERROR' $WORK_CTL/$CTL_LOG
      status=$?
      if ((status == 0)); then
        echo "** PROBLEM WITH CTL RUN. CHECK LOG FILE."
        save_ctl=1
        failed=1
      fi

# did the i/j to lat/lon transform fail for the test?
      grep -Eq 'BAD|ERROR' $WORK_TEST/$TEST_LOG
      status=$?
      if ((status == 0)); then
        echo "** PROBLEM WITH TEST RUN. CHECK LOG FILE."
        save_test=1
        failed=1
      fi

# if any step of regression test fails, save data in a subdirectory.

      if ((save_ctl == 1));then
        FAILED_DIR=$WORK_CTL/failed.${bytesize}byte.grid${grids}.${routine}
        mkdir -p $FAILED_DIR
        if [ -s $WORK_CTL/$CTL_LOG ]; then
          mv $WORK_CTL/$CTL_LOG $FAILED_DIR
        fi
        if [ -s $WORK_CTL/grid${grids}.bin ];then
          mv $WORK_CTL/grid${grids}.bin $FAILED_DIR
        fi
        if [ -s $WORK_CTL/core* ];then
          mv $WORK_CTL/core*  $FAILED_DIR
        fi 
      fi

      if ((save_test == 1));then
        FAILED_DIR=$WORK_TEST/failed.${bytesize}byte.grid${grids}.${routine}
        mkdir -p $FAILED_DIR
        if [ -s $WORK_TEST/$TEST_LOG ];then
          mv $WORK_TEST/$TEST_LOG $FAILED_DIR
        fi
        if [ -s $WORK_TEST/grid${grids}.bin ];then
          mv $WORK_TEST/grid${grids}.bin $FAILED_DIR
        fi
        if [ -s $WORK_TEST/core* ];then
          mv $WORK_TEST/core*  $FAILED_DIR
        fi
      fi

      rm -f $WORK_CTL/grid${grids}.bin $WORK_TEST/grid${grids}.bin
      rm -f $WORK_CTL/$CTL_LOG $WORK_TEST/$TEST_LOG
      rm -f $WORK_CTL/core*  $WORK_TEST/core*

    done
  done
done

if ((failed == 0)); then
  echo
  echo "<<< GDSWIZ/WZD REGRESSION TEST PASSED. >>>"
  echo
else
  echo
  echo "<<< GDSWIZ/WZD REGRESSION TEST FAILED. >>>"
  echo
fi

exit 0
