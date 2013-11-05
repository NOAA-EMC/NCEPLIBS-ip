#!/bin/ksh

#--------------------------------------------------------------
# Run regression test for iplib routine ipmerge2.
#
# The routine is invoked by a simple Fortran program.
# The program is compiled with all three byte versions
# of the 'control' and 'test' ip library.
# 
# The three byte versions of the library are:
#  > 4 byte integer/4 byte float  ($bytesize=4)
#  > 8 byte integer/8 byte float  ($bytesize=8)
#  > 8 byte float/4 byte integer  ($bytesize=d)
#
# Output from the program is written to a text log file.
#
# The log file naming convention is:
#    ctl_${bytesize}byte.log 
#    test_${bytesize}byte.log
#
# The 'control' and 'test' libraries must produce identical
# output, or the regression test is considered 'failed'. 
# If a failure happens, the log files are stored in the work 
# directory with a ".failed" extension.
# 
# This script is run by the Runall.${machine}.ksh driver 
# script located in /reg_tests.
#--------------------------------------------------------------

#set -x

echo
echo BEGIN IPMERGE2 REGRESSION TEST
echo

REG_DIR=${REG_DIR:-../..}

WORK_DIR=${WORK_DIR:-/stmp/$LOGNAME/regression}

# where the executables are located.
EXEC_DIR=$REG_DIR/ipmerge2/exec

WORK=${WORK_DIR}/ipmerge2
rm -fr $WORK
mkdir -p $WORK
WORK_CTL=${WORK}/ctl
mkdir -p $WORK_CTL
cp $EXEC_DIR/ipmerge2_ctl_*.exe  $WORK_CTL
WORK_TEST=${WORK}/test
mkdir -p $WORK_TEST
cp $EXEC_DIR/ipmerge2_test_*.exe $WORK_TEST

failed=0

for bytesize in "4" "8" "d"  # test all versions of library
do

  echo TEST ${bytesize}-BYTE FLOAT VERSION OF ROUTINE IPMERGE2

  cd $WORK_CTL
  CTL_LOG=ctl.${bytesize}byte.log
  ipmerge2_ctl_${bytesize}.exe  > $CTL_LOG
  status=$?
  if ((status != 0)); then
    echo CONTROL RUN FAILED.
    failed=1
    continue
  fi

  cd $WORK_TEST
  TEST_LOG=test.${bytesize}byte.log
  ipmerge2_test_${bytesize}.exe > $TEST_LOG
  status=$?
  if ((status != 0)); then
    echo TEST RUN FAILED.
    failed=1
    continue
  fi

  cmp $WORK_CTL/$CTL_LOG $WORK_TEST/$TEST_LOG
  status=$?
  if ((status != 0)); then
    echo LOG FILES NOT BIT IDENTIAL. REGRESSION TEST FAILED.
    echo CHECK LOG FILES SAVED IN WORK DIRECTORY.
    mv $WORK_CTL/$CTL_LOG  $WORK_CTL/${CTL_LOG}.failed
    mv $WORK_TEST/$TEST_LOG $WORK_TEST/${TEST_LOG}.failed
    failed=1
  fi

done

if ((failed == 0));then
  echo
  echo "<<< IPMERGE2 REGRESSION TEST PASSED. >>>"
  echo
else
  echo
  echo "<<< IPMERGE2 REGRESSION TEST FAILED. >>>"
  echo
fi

exit 0
