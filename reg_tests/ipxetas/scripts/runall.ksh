#!/bin/ksh

#------------------------------------------------------------------------
# Run a Fortran program to test iplib routine ipxetas as follows: 
#
# 1) Convert land/sea mask on a rotated lat/lon "e"-staggered "h" point 
#    grid to a rotated lat/lon unstaggerd "full" grid. 
# 2) Convert u-component wind on a rotated lat/lon "e"-staggered "v" point 
#    grid to a rotated lat/lon unstaggerd "full" grid. 
# 3) Convert land/sea mask on a rotated lat/lon unstaggered "full" grid
#    to a rotated lat/lon "e"-staggered "h" point grid. 
# 4) Convert land/sea mask on a rotated lat/lon unstaggered "full" grid
#    to a rotated lat/lon "e"-staggered "v" point grid. 
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
# The input files are in the ../data directory.  They are grib 2 format.
#
# Output files from the program are in grib 2 format.  If the files from
# the 'test' and 'control' ip libraries are not bit identical the 
# the test has "failed".  And the output grib 2 files are saved in the
# working directory and named as:
#
#   output.${cases}.${bytesize}byte.failed.grb2
#
# where ${cases}="h2full" or "v2full" or "full2h" or "full2v"
#
# Also, the program log files from failed tests are saved in the
# working directories and named:
#
#   ${cases}.${bytesize}byte.log.failed
#
# This script is run by the Runall.${machine}.ksh driver script located
# in /reg_tests.
#------------------------------------------------------------------------

#set -x

echo
echo BEGIN REGRESSION TEST FOR IPXETAS
echo

WORK_DIR=${WORK_DIR:-/stmpp1/$LOGNAME/regression}

REG_DIR=${REG_DIR:-$(pwd)/../..}

EXEC_DIR=${REG_DIR}/ipxetas/exec

DATA_DIR=${REG_DIR}/ipxetas/data

WORK=${WORK_DIR}/ipxetas
rm -fr $WORK
mkdir -p $WORK

WORK_CTL=${WORK}/ctl
mkdir -p $WORK_CTL
cp $EXEC_DIR/ipxetas_ctl_*.exe  $WORK_CTL

WORK_TEST=${WORK}/test
mkdir -p $WORK_TEST
cp $EXEC_DIR/ipxetas_test_*.exe $WORK_TEST

reg_test_failed=0

for bytesize in "4" "8" "d"  # the three byte versions of the library
do

  for cases in "full2h" "full2v" "h2full" "v2full"
  do

    echo TEST ${bytesize}-BYTE VERSION OF IPXETAS FOR CASE: $cases

    case $cases in
    "full2h")
      INPUT_DATA=$DATA_DIR/full_slmask.grb2
      IDIR=-1 ;;
    "full2v")
      INPUT_DATA=$DATA_DIR/full_slmask.grb2
      IDIR=-2 ;;
    "h2full")
      INPUT_DATA=$DATA_DIR/egrid_hpnt_slmask.grb2
      IDIR=0 ;;
    "v2full")
      INPUT_DATA=$DATA_DIR/egrid_vpnt_uwind.grb2
      IDIR=0 ;;
    esac

    cp $INPUT_DATA $WORK_CTL/fort.9
    cp $INPUT_DATA $WORK_TEST/fort.9

    ctl_failed=0
    test_failed=0

    cd $WORK_CTL
    CTL_LOG=${cases}.${bytesize}byte.log
    ipxetas_ctl_${bytesize}.exe $IDIR > $CTL_LOG
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
    TEST_LOG=${cases}.${bytesize}byte.log
    ipxetas_test_${bytesize}.exe $IDIR > $TEST_LOG
    status=$?
    if ((status != 0)); then
      echo "** TEST RUN FAILED."
      test_failed=1
      reg_test_failed=1
      if [ -s $WORK_TEST/$TEST_LOG ];then
        mv $WORK_TEST/$TEST_LOG  $WORK_TEST/${TEST_LOG}.failed
      fi
    fi

# If test and control executables ran to completion, check output
# grib 2 files for bit identicalness.

    if ((ctl_failed == 0 && test_failed == 0));then

      cmp $WORK_CTL/output.grb2 $WORK_TEST/output.grb2
      status=$?
      if ((status != 0)); then
        echo "** OUTPUT GRIB 2 FILES NOT BIT IDENTICAL. REGRESSION TEST FAILED."
        mv $WORK_CTL/output.grb2 $WORK_CTL/output.${cases}.${bytesize}byte.failed.grb2
        mv $WORK_TEST/output.grb2 $WORK_TEST/output.${cases}.${bytesize}byte.failed.grb2
        mv $WORK_CTL/$CTL_LOG  $WORK_CTL/${CTL_LOG}.failed
        mv $WORK_TEST/$TEST_LOG  $WORK_TEST/${TEST_LOG}.failed
        reg_test_failed=1
      fi

    fi 

    rm -f $WORK_CTL/$CTL_LOG  $WORK_TEST/$TEST_LOG
    rm -f $WORK_CTL/output.grb2 $WORK_TEST/output.grb2

  done

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
