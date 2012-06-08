#!/bin/ksh

#--------------------------------------------------------------
# Test iplib routine ipmerge2, which merges two bitmap 
# fields.
# 
# Output from the 'control' and 'test' is placed in its
# own text log file.  If the log files are not identical,
# the test is considered 'failed'.  If a test fails, the
# log file is stored in the work directory with a ".failed"
# extension.
# 
# All three versions of the library are tested:
#  > 4 byte integer/4 byte float  (libip_4.a)
#  > 8 byte integer/8 byte float  (libip_8.a)
#  > 8 byte float/4 byte integer  (libip_d.a)
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
cp $EXEC_DIR/ctl/*exe $WORK_CTL
WORK_TEST=${WORK}/test
mkdir -p $WORK_TEST
cp $EXEC_DIR/test/*exe $WORK_TEST

failed=0

for bytesize in "4" "8" "d"  # test all versions of library
do
  echo TEST ${bytesize}-BYTE FLOAT VERSION OF ROUTINE IPMERGE2
  cd $WORK_CTL
  CTL_LOG=ctl.${bytesize}byte.log
  ipmerge2_ctl_${bytesize}.exe  > $CTL_LOG
  cd $WORK_TEST
  TEST_LOG=test.${bytesize}byte.log
  ipmerge2_test_${bytesize}.exe > $TEST_LOG
  cmp $WORK_CTL/$CTL_LOG $WORK_TEST/$TEST_LOG
  status=$?
  if ((status != 0))
  then
    echo LOG FILES NOT BIT IDENTIAL. TEST FAILED.
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
