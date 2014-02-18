#!/bin/ksh

#------------------------------------------------------------------------
# Test iplib routine ipxetas as follows: run a Fortran program
# to read an input file of vegetation greenness on a 'full' 12km eta
# grid, then call routine ipxetas to:
#
# 1) create a staggered mass grid from the full grid.
# 2) create a staggered velocity grid from the full grid.
# 3) create a full grid from the staggered mass grid created by step (1)
# 4) create a full grid from the staggered vel grid created by step (2)
#
# The Fortran program is compiled with all three byte versions
# of the 'control' and 'test' ip library.
#
# The three byte versions of the library are:
#  > 4 byte integer/4 byte float  ($bytesize=4)
#  > 8 byte integer/8 byte float  ($bytesize=8)
#  > 8 byte float/4 byte integer  ($bytesize=d)
#
# The 'control' and 'test' executables run in their own working directories:
# $WORK_CTL and $WORK_TEST.
#
# The input greenness data is: ../data/green.202.grb.  It is in
# grib 1 format.
#
# Output from steps (1) and (2) is in the binary file named "staggered.bin"
# Output from steps (3) and (4) is in the binary file named "full.bin"
#
# If the binary files from the 'test' and 'control' ip libraries are not bit
# identical the test has "failed".  And the binary files are renamed as:
#
#   staggered.${bytesize}byte.bin.failed
#   full.${bytesize}byte.bin.failed
#
# The routine also computes a kgds array used by the w3nco library.
# This array is output to a text log file.  If the 'test' and 'control' 
# log files are not identical, the regression test is considered "failed".
# And the log files are renamed as:
#
#   ctl.${bytesize}byte.log.failed
#   test.${bytesize}byte.log.failed
#
# This script is run by the Runall.${machine}.ksh driver script located
# in /reg_tests.
#------------------------------------------------------------------------

#set -x

echo
echo BEGIN REGRESSION TEST FOR IPXETAS
echo

WORK_DIR=${WORK_DIR:-/stmpp1/$LOGNAME/regression}

REG_DIR=${REG_DIR:-../..}

EXEC_DIR=${REG_DIR}/ipxetas/exec
INPUT_DATA=$REG_DIR/ipxetas/data/green.202.grb

WORK=${WORK_DIR}/ipxetas
rm -fr $WORK
mkdir -p $WORK
WORK_CTL=${WORK}/ctl
mkdir -p $WORK_CTL
cp $EXEC_DIR/ipxetas_ctl_*.exe  $WORK_CTL
cp $INPUT_DATA $WORK_CTL/fort.9
WORK_TEST=${WORK}/test
mkdir -p $WORK_TEST
cp $EXEC_DIR/ipxetas_test_*.exe $WORK_TEST
cp $INPUT_DATA $WORK_TEST/fort.9

reg_test_failed=0

for bytesize in "4" "8" "d"  # the three byte versions of the library
do

  echo TEST ${bytesize}-BYTE VERSION OF IPXETAS

  ctl_failed=0
  test_failed=0

  cd $WORK_CTL
  CTL_LOG=ctl.${bytesize}byte.log
  ipxetas_ctl_${bytesize}.exe  > $CTL_LOG
  status=$?
  if ((status != 0)); then
    echo "** CONTROL RUN FAILED."
    ctl_failed=1
    reg_test_failed=1
    if [ -s $WORK_CTL/$CTL_LOG ];then
      mv $WORK_CTL/$CTL_LOG  $WORK_CTL/${CTL_LOG}.failed
    fi
  fi

  cd $WORK_TEST
  TEST_LOG=test.${bytesize}byte.log
  ipxetas_test_${bytesize}.exe > $TEST_LOG
  status=$?
  if ((status != 0)); then
    echo "** TEST RUN FAILED."
    test_failed=1
    reg_test_failed=1
    if [ -s $WORK_TEST/$TEST_LOG ];then
      mv $WORK_TEST/$TEST_LOG  $WORK_TEST/${TEST_LOG}.failed
    fi
  fi

# if test and control executables ran to completion, check binary 
# and log files for bit identicalness.

  if ((ctl_failed == 0 && test_failed == 0));then

    save_log=0

    cmp $WORK_CTL/staggered.bin $WORK_TEST/staggered.bin
    status=$?
    if ((status != 0)); then
      echo "** STAGGERED GRID BINARY FILES NOT BIT IDENTICAL. REGRESSION TEST FAILED."
      mv $WORK_CTL/staggered.bin $WORK_CTL/staggered.${bytesize}byte.bin.failed
      mv $WORK_TEST/staggered.bin $WORK_TEST/staggered.${bytesize}byte.bin.failed
      reg_test_failed=1
      save_log=1
    fi

    cmp $WORK_CTL/full.bin $WORK_TEST/full.bin
    status=$?
    if ((status != 0)); then
      echo "** FULL GRID BINARY FILES NOT BIT IDENTICAL. REGRESSION TEST FAILED."
      mv $WORK_CTL/full.bin $WORK_CTL/full.${bytesize}byte.bin.failed
      mv $WORK_TEST/full.bin $WORK_TEST/full.${bytesize}byte.bin.failed
      reg_test_failed=1
      save_log=1
    fi

    cmp $WORK_CTL/$CTL_LOG  $WORK_TEST/$TEST_LOG
    status=$?
    if ((status != 0)); then
      echo "** KGDS ARRAYS NOT BIT IDENTICAL. REGRESSION TEST FAILED."
      reg_test_failed=1
      save_log=1
    fi

    if ((save_log == 1));then
      mv $WORK_CTL/$CTL_LOG  $WORK_CTL/${CTL_LOG}.failed
      mv $WORK_TEST/$TEST_LOG  $WORK_TEST/${TEST_LOG}.failed
    fi

  fi 

  rm -f $WORK_CTL/$CTL_LOG  $WORK_TEST/$TEST_LOG
  rm -f $WORK_CTL/staggered.bin $WORK_TEST/staggered.bin
  rm -f $WORK_CTL/full.bin $WORK_TEST/full.bin

done

if ((reg_test_failed == 0));then
  echo
  echo "<<< IPXETAS REGRESSION TEST PASSED. >>>"
  echo
else
  echo
  echo "<<< IPXETAS REGRESSION TEST FAILED. >>>"
  echo
fi


exit 0
