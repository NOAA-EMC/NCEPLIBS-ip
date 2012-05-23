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

for bytesize in "4" "8" "d"  # all three byte versions of the library.
do
  save_ctl_log=0
  save_test_log=0
  echo TEST ${bytesize}-BYTE FLOAT VERSION OF ROUTINE MAKGDS
  cd $WORK_CTL
  makgds_ctl_${bytesize}.exe  > ctl.log
  grep -Eq 'BAD|ERROR' $WORK_CTL/ctl.log
  status=$?
  if ((status == 0)); then
    echo PROBLEM WITH CTL RUN. CHECK LOG FILE.
    save_ctl_log=1
  fi
  cd $WORK_TEST
  makgds_test_${bytesize}.exe > test.log
  grep -Eq 'BAD|ERROR' $WORK_TEST/test.log
  status=$?
  if ((status == 0)); then
    echo PROBLEM WITH TEST RUN. CHECK LOG FILE.
    save_test_log=1
  fi
  cmp $WORK_CTL/ctl.log $WORK_TEST/test.log
  status=$?
  if ((status != 0))
  then
    echo LOG FILES NOT BIT IDENTIAL. TEST FAILED.
    save_ctl_log=1
    save_test_log=1
  fi

  if ((save_ctl_log == 1)); then
    mv $WORK_CTL/ctl.log $WORK_CTL/ctl.${bytesize}byte.log.failed
  fi

  if ((save_test_log == 1)); then
    mv $WORK_TEST/test.log $WORK_TEST/test.${bytesize}byte.log.failed
  fi

done

echo
echo MAKGDS REGRESSION TEST COMPLETED.
echo

exit 0
