#!/bin/ksh

#-------------------------------------------------------------------
# Regression test for iplib routine makgds.
#
# Routine is tested as follows:
#  1) create gds and kgds arrays for ncep grid 3.  arrays hold grid
#     description information used by w3 grib library.
#  2) make kgds array for grid 3 from gds array
#  3) make gds array for grid 3 from kgds array
#
# Output from the control and test executables is placed in a
# ascii log file.  If the log files are not bit identical, the
# test fails. If a test fails, the log file is saved in 
# WORK_DIR with a ".failed" extension.
#-------------------------------------------------------------------

#set -x

echo
echo BEGIN MAKGDS REGRESSION TEST
echo
 
WORK_DIR=${WORK_DIR:-/stmp/$LOGNAME/regression}

REG_DIR=${REG_DIR:-../..}

EXEC_DIR=$REG_DIR/makgds/exec

WORK=$WORK_DIR/makgds
rm -fr $WORK
mkdir -p $WORK
WORK_CTL=${WORK}/ctl
mkdir -p $WORK_CTL
cp $EXEC_DIR/ctl/*exe $WORK_CTL
WORK_TEST=${WORK}/test
mkdir -p $WORK_TEST
cp $EXEC_DIR/test/*exe $WORK_TEST

reg_test_failed=0

for bytesize in "4" "8" "d"  # all three byte versions of the library.
do

  ctl_failed=0
  test_failed=0

  echo TEST ${bytesize}-BYTE FLOAT VERSION OF ROUTINE MAKGDS

  cd $WORK_CTL
  makgds_ctl_${bytesize}.exe  > ctl.log
  status=$?
  if ((status != 0)); then
    echo "** PROBLEM WITH CTL RUN. CHECK LOG FILE."
    ctl_failed=1
    reg_test_failed=1
  fi

  cd $WORK_TEST
  makgds_test_${bytesize}.exe > test.log
  status=$?
  if ((status != 0)); then
    echo "** PROBLEM WITH TEST RUN. CHECK LOG FILE."
    test_failed=1
    reg_test_failed=1
  fi

  if ((ctl_failed == 0 && test_failed == 0));then
    cmp $WORK_CTL/ctl.log $WORK_TEST/test.log
    status=$?
    if ((status != 0))
    then
      echo "** LOG FILES NOT BIT IDENTICAL. REGRESSION TEST FAILED."
      reg_test_failed=1
      ctl_failed=1
      test_failed=1
    fi
  fi

  if ((ctl_failed == 1));then
    if [ -s $WORK_CTL/ctl.log ];then
      mv $WORK_CTL/ctl.log $WORK_CTL/ctl.${bytesize}byte.log.failed
    fi
  fi

  if ((test_failed == 1));then
    if [ -s $WORK_TEST/test.log ];then
      mv $WORK_TEST/test.log $WORK_TEST/test.${bytesize}byte.log.failed
    fi
  fi

  rm -f $WORK_TEST/test.log  $WORK_CTL/ctl.log

done

if ((reg_test_failed == 0)); then
  echo
  echo "<<< MAKGDS REGRESSION TEST PASSED. >>>"
  echo
else
  echo
  echo "<<< MAKGDS REGRESSION TEST FAILED. >>>"
  echo
fi

exit 0
