#!/bin/ksh

#--------------------------------------------------------------
# Test iplib routines gcdist and movect, which 
# compute fields associated with great circle routes.
#
# Output from the 'control' and 'test' is placed in its
# own text log file.  If the log files are not identical,
# the test is considered 'failed'.
#
# All three versions of the library are tested:
#  > 4 byte integer/4 byte float  (libip_4.a)
#  > 8 byte integer/8 byte float  (libip_8.a)
#  > 8 byte float/4 byte integer  (libip_d.a)
#--------------------------------------------------------------

#set -x

echo
echo "BEGIN GCDIST/MOVECT REGRESSION TEST"
echo

REG_DIR=${REG_DIR:-../..}

WORK_DIR=${WORK_DIR:-/stmp/$LOGNAME/regression}

# where the executables are located
EXEC_DIR=$REG_DIR/gcdist/exec

WORK=${WORK_DIR}/gcdist
rm -fr $WORK
mkdir -p $WORK
WORK_CTL=${WORK}/ctl
mkdir -p $WORK_CTL
cp $EXEC_DIR/ctl/*exe $WORK_CTL
WORK_TEST=${WORK}/test
mkdir -p $WORK_TEST
cp $EXEC_DIR/test/*exe $WORK_TEST

failed=0

for bytesize in "4" "8" "d"  # the version of the library
do

  echo "TEST ${bytesize}-BYTE FLOAT VERSION OF ROUTINES GCDIST/MOVECT"

  cd $WORK_CTL
  CTL_LOG=ctl.${bytesize}byte.log
  gcdist_ctl_${bytesize}.exe  > $CTL_LOG
  status=$?
  if ((status != 0))
  then
    echo CONTROL FAILED.
    failed=1
  fi

  cd $WORK_TEST
  TEST_LOG=test.${bytesize}byte.log
  gcdist_test_${bytesize}.exe > $TEST_LOG
  status=$?
  if ((status != 0))
  then
    echo TEST FAILED.
    failed=1
  fi

  cmp $WORK_CTL/$CTL_LOG $WORK_TEST/$TEST_LOG
  status=$?
  if ((status != 0))
  then
    echo "LOG FILES NOT BIT IDENTIAL. TEST FAILED."
    echo "CHECK LOG FILES SAVED IN WORK DIRECTORY."
    mv $WORK_CTL/$CTL_LOG $WORK_CTL/${CTL_LOG}.failed
    mv $WORK_TEST/$TEST_LOG $WORK_TEST/${TEST_LOG}.failed
    failed=1
  fi

done

if ((failed == 0))
then
  echo
  echo "<<< GCDIST/MOVECT REGRESSION TEST PASSED. >>>"
  echo
else
  echo
  echo "<<< GCDIST/MOVECT REGRESSION TEST FAILED. >>>"
  echo
fi

exit 0
