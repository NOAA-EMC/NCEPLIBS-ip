#!/bin/ksh

#--------------------------------------------------------------
# Test iplib routines ipsector and ipspaste, which
# creates a subset of a larger two-dimensional field,
# and routine ipspaste, which does the opposite.
#
# Output from each step of the 'control' and 'test' is placed
# in its own binary file (named ipsector.bin and ipspaste.bin).  
# If the files are not identical, the test is considered 'failed'.
# And the binary files are saved with a ".failed" extension.
#
# All three versions of the library are tested:
#  > 4 byte integer/4 byte float  (libip_4.a)
#  > 8 byte integer/8 byte float  (libip_8.a)
#  > 8 byte float/4 byte integer  (libip_d.a)
#--------------------------------------------------------------

#set -x

echo
echo BEGIN IPSECTOR/IPSPASTE REGRESSION TEST
echo

WORK_DIR=${WORK_DIR:-/stmp/$LOGNAME/regression}

REG_DIR=${REG_DIR:-../..}

EXEC_DIR=$REG_DIR/ipsector/exec
INPUT_DATA=$REG_DIR/ipsector/data/global_tg3clim.1x1.grb

WORK=${WORK_DIR}/ipsector
rm -fr $WORK
mkdir -p $WORK
WORK_CTL=${WORK}/ctl
mkdir -p $WORK_CTL
cp $EXEC_DIR/ctl/*exe $WORK_CTL
cp $INPUT_DATA $WORK_CTL/fort.9
WORK_TEST=${WORK}/test
mkdir -p $WORK_TEST
cp $EXEC_DIR/test/*exe $WORK_TEST
cp $INPUT_DATA $WORK_TEST/fort.9

for bytesize in "4" "8" "d"  # the byte version of the ip library
do

  echo TEST ${bytesize}-BYTE VERSION OF IPSECTOR/IPSPASTE
  cd $WORK_CTL
  CTL_LOG=ctl.${bytesize}byte.log
  ipsector_ctl_${bytesize}.exe  > $CTL_LOG
  cd $WORK_TEST
  TEST_LOG=test.${bytesize}byte.log
  ipsector_test_${bytesize}.exe > $TEST_LOG

  save_ctl_log=0
  save_test_log=0

  cmp $WORK_CTL/ipsector.bin $WORK_TEST/ipsector.bin
  status=$?
  if ((status != 0))
  then
    echo IPSECTOR BINARY FILES NOT BIT IDENTICAL. TEST FAILED.
    mv $WORK_CTL/ipsector.bin $WORK_CTL/ipsector.${bytesize}byte.bin.failed
    mv $WORK_TEST/ipsector.bin $WORK_TEST/ipsector.${bytesize}byte.bin.failed
    save_ctl_log=1
    save_test_log=1
  fi

  cmp $WORK_CTL/ipspaste.bin $WORK_TEST/ipspaste.bin
  status=$?
  if ((status != 0))
  then
    echo IPSPASTE BINARY FILES NOT BIT IDENTICAL. TEST FAILED.
    mv $WORK_CTL/ipspaste.bin $WORK_CTL/ipspaste.${bytesize}byte.bin.failed
    mv $WORK_TEST/ipspaste.bin $WORK_TEST/ipspaste.${bytesize}byte.bin.failed
    save_ctl_log=1
    save_test_log=1
  fi

  grep -Eq 'BAD|ERROR' $WORK_CTL/$CTL_LOG
  status=$?
  if ((status == 0)); then
    echo PROBLEM WITH CTL RUN. CHECK LOG FILE.
    save_ctl_log=1
  fi

  grep -Eq 'BAD|ERROR' $WORK_TEST/$TEST_LOG
  status=$?
  if ((status == 0)); then
    echo PROBLEM WITH TEST RUN. CHECK LOG FILE.
    save_test_log=1
  fi

  if ((save_ctl_log == 1)); then
    mv $WORK_CTL/$CTL_LOG $WORK_CTL/${CTL_LOG}.failed
  fi

  if ((save_test_log == 1)); then
    mv $WORK_TEST/$TEST_LOG $WORK_TEST/${TEST_LOG}.failed
  fi

  rm -f $WORK_CTL/ipsector.bin $WORK_TEST/ipsector.bin
  rm -f $WORK_CTL/ipspaste.bin $WORK_TEST/ipspaste.bin
done

echo
echo IPSECTOR/IPSPASTE REGRESSION TEST COMPLETED.
echo

exit 0
