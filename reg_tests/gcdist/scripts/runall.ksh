#!/bin/ksh

#--------------------------------------------------------------
# Test iplib routines gcdist and movect, which 
# compute fields associated with great circle routes.
#
# Output from the 'control' and 'test' is placed in its
# own text log file.  If the log files are not identical,
# the regression test is considered 'failed'.
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

reg_test_failed=0

for bytesize in "4" "8" "d"  # the version of the library
do

  echo "TEST ${bytesize}-BYTE FLOAT VERSION OF ROUTINES GCDIST/MOVECT"

  ctl_failed=0
  test_failed=0

  cd $WORK_CTL
  CTL_LOG=ctl.${bytesize}byte.log
  gcdist_ctl_${bytesize}.exe  > $CTL_LOG
  status=$?
  if ((status != 0));then
    echo "** CONTROL RUN FAILED."
    reg_test_failed=1
    ctl_failed=1
  fi

  cd $WORK_TEST
  TEST_LOG=test.${bytesize}byte.log
  gcdist_test_${bytesize}.exe > $TEST_LOG
  status=$?
  if ((status != 0));then
    echo "** TEST RUN FAILED."
    reg_test_failed=1
    test_failed=1
  fi

  if ((ctl_failed == 0 && test_failed == 0));then
    cmp $WORK_CTL/$CTL_LOG $WORK_TEST/$TEST_LOG
    status=$?
    if ((status != 0));then
      echo "** LOG FILES NOT BIT IDENTIAL. REGRESSION TEST FAILED."
      echo "** CHECK LOG FILES SAVED IN WORK DIRECTORY."
      reg_test_failed=1
      ctl_failed=1
      test_failed=1
    fi
  fi

  if ((ctl_failed == 1));then
    if [ -s $WORK_CTL/$CTL_LOG ];then
      mv $WORK_CTL/$CTL_LOG $WORK_CTL/${CTL_LOG}.failed
    fi
  fi

  if ((test_failed == 1));then
    if [ -s $WORK_TEST/$TEST_LOG ];then
      mv $WORK_TEST/$TEST_LOG $WORK_TEST/${TEST_LOG}.failed
    fi
  fi

done

if ((reg_test_failed == 0))
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
