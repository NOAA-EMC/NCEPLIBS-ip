#!/bin/ksh

#--------------------------------------------------------------
# Test iplib routines ipsector and ipspaste, which
# creates a subset of a larger two-dimensional field,
# and routine ipspaste, which does the opposite.
#
# Read in a global dataset of substrate
# temperature, then call ipsector to create a subset of 
# the original grid.  Then call ipspaste to 'paste'
# the subset data (created by call to ipsector) back
# to the original grid.  The data returned from ipsector
# should match the original data.
#
# Regression test passes if...
#  - The 'control' and 'test' runs to completion
#    with no errors. 
#  - The output from the call to ipsector is stored in
#    file ipsector.bin. The 'control' and 'test' 
#    files must be bit identical.
#  - The output from the call to ipspaste is stored in
#    file ipspaste.bin. The 'control' and 'test' 
#    files must be bit identical.
#  - The output from ipspaste must be bit identical
#    to the input data (stored in input.bin).
#    Check performed for 'test' files only.
#  - The kgds array values from the 'control' and
#    'test' are stored to log files.  The log
#    files must be identical.
#
# If anything fails, the log and binary data are stored
# in the work directory with a 'failed' extension.
#
# These ip routines have separate logic for differing
# Grib 1 scanning modes.  So there are two versions
# of the input data - one with a scan mode of '0'
# and one with a scan mode of '32'.
#
# All three byte versions of the library are tested:
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

WORK=${WORK_DIR}/ipsector
rm -fr $WORK
mkdir -p $WORK
WORK_CTL=${WORK}/ctl
mkdir -p $WORK_CTL
cp $EXEC_DIR/ctl/*exe $WORK_CTL
WORK_TEST=${WORK}/test
mkdir -p $WORK_TEST
cp $EXEC_DIR/test/*exe $WORK_TEST

INPUT_DATA=$REG_DIR/ipsector/data/global_tg3clim.1x1.scan.0.grb
cp $INPUT_DATA $WORK_CTL
cp $INPUT_DATA $WORK_TEST
INPUT_DATA=$REG_DIR/ipsector/data/global_tg3clim.1x1.scan.32.grb
cp $INPUT_DATA $WORK_CTL
cp $INPUT_DATA $WORK_TEST

for scan in "0" "32"   # grib 1 scanning mode 
do

  ln -fs $WORK_CTL/global_tg3clim.1x1.scan.${scan}.grb   $WORK_CTL/fort.9
  ln -fs $WORK_TEST/global_tg3clim.1x1.scan.${scan}.grb  $WORK_TEST/fort.9

  for bytesize in "4" "8" "d"  # the byte version of the ip library
  do

    echo "TEST ${bytesize}-BYTE VERSION OF IPSECTOR/IPSPASTE FOR SCAN MODE ${scan}"

    save_ctl=0
    save_test=0

    cd $WORK_CTL
    CTL_LOG=ctl.${bytesize}byte.log
    ipsector_ctl_${bytesize}.exe  > $CTL_LOG
    status=$?
    if ((status != 0))
    then
      echo "CONTROL RUN FAILED."
      save_ctl=1
    fi

    cd $WORK_TEST
    TEST_LOG=test.${bytesize}byte.log
    ipsector_test_${bytesize}.exe > $TEST_LOG
    status=$?
    if ((status != 0))
    then
      echo "TEST RUN FAILED."
      save_test=1
    fi

    cmp $WORK_CTL/ipsector.bin $WORK_TEST/ipsector.bin
    status=$?
    if ((status != 0))
    then
      echo "IPSECTOR BINARY FILES NOT BIT IDENTICAL. TEST FAILED."
      save_ctl=1
      save_test=1
    fi

    cmp $WORK_CTL/ipspaste.bin $WORK_TEST/ipspaste.bin
    status=$?
    if ((status != 0))
    then
      echo "IPSPASTE BINARY FILES NOT BIT IDENTICAL. TEST FAILED."
      save_ctl=1
      save_test=1
    fi

    cmp $WORK_TEST/input.bin $WORK_TEST/ipspaste.bin
    status=$?
    if ((status != 0))
    then
      echo "IPSPASTE BINARY FILE AND INPUT FILE NOT BIT IDENTICAL. TEST FAILED."
      save_test=1
    fi

    cmp $WORK_CTL/$CTL_LOG $WORK_TEST/$TEST_LOG
    status=$?
    if ((status != 0))
    then
      echo "KGDS ARRAY NOT IDENTICAL. TEST FAILED."
      save_ctl=1
      save_test=1
    fi

    grep -Eq 'BAD|ERROR' $WORK_CTL/$CTL_LOG
    status=$?
    if ((status == 0)); then
      echo "PROBLEM WITH CTL RUN. TEST FAILED."
      save_ctl=1
    fi

    grep -Eq 'BAD|ERROR' $WORK_TEST/$TEST_LOG
    status=$?
    if ((status == 0)); then
      echo "PROBLEM WITH TEST RUN. TEST FAILED."
      save_test=1
    fi

    if ((save_ctl == 1)); then
      if [ -s $WORK_CTL/$CTL_LOG ];then
        mv $WORK_CTL/$CTL_LOG $WORK_CTL/${CTL_LOG}.${bytesize}byte.${scan}scan.failed
      fi
      if [ -s $WORK_CTL/ipspaste.bin ];then
        mv $WORK_CTL/ipspaste.bin $WORK_CTL/ipspaste.${bytesize}byte.${scan}scan.failed.bin
      fi
      if [ -s $WORK_CTL/ipsector.bin ];then
        mv $WORK_CTL/ipsector.bin $WORK_CTL/ipsector.${bytesize}byte.${scan}scan.failed.bin
      fi
    fi

    if ((save_test == 1)); then
      if [ -s $WORK_TEST/$TEST_LOG ];then
        mv $WORK_TEST/$TEST_LOG $WORK_TEST/${TEST_LOG}.${bytesize}byte.${scan}scan.failed
      fi 
      if [ -s $WORK_TEST/ipspaste.bin ];then
        mv $WORK_TEST/ipspaste.bin $WORK_TEST/ipspaste.${bytesize}byte.${scan}scan.failed.bin
      fi
      if [ -s $WORK_TEST/ipsector.bin ];then
        mv $WORK_TEST/ipsector.bin $WORK_TEST/ipsector.${bytesize}byte.${scan}scan.failed.bin
      fi
    fi

    rm -f $WORK_CTL/ipsector.bin $WORK_TEST/ipsector.bin
    rm -f $WORK_CTL/ipspaste.bin $WORK_TEST/ipspaste.bin

  done  # library byte size loop

  rm -f $WORK_CTL/fort.9 $WORK_TEST/fort.9
  echo

done # scan mode loop

echo
echo "IPSECTOR/IPSPASTE REGRESSION TEST COMPLETED."
echo

exit 0
