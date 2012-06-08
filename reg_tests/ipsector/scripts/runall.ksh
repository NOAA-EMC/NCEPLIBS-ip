#!/bin/ksh

#--------------------------------------------------------------
# Test iplib routine ipsector, which creates
# a subsector of a larger two-dimensional field,
# and routine ipspaste, which does the opposite.
#
# Read in a global dataset of substrate temperature,
# then call ipsector to create a subsector of 
# the original grid.  Then call ipspaste to 'paste'
# the subsectored data (created by call to ipsector) back
# to the original grid.  The data returned from ipsector
# should match the original data.
#
# Three sets of calls to ipsector/ipspaste are made:
# - for a North America subsector
# - for a non-overlapping subsector
# - for an overlapping subsector
#
# Regression test passes if...
#  - The 'control' and 'test' runs to completion
#    with no errors. 
#  - The output from the three calls to ipsector are stored in
#    files "ipsector.namer.bin", "ipsector.no.sect.bin" and
#    "ipsector.ovlp.sect.bin". The 'control' and 'test'
#    files must be bit identical.
#  - The output from the three calls to ipspaste are stored in
#    files "ipspaste.namer.bin", "ipspaste.no.sect.bin" and
#    "ipspaste.ovlp.sect.bin". These files must
#    be bit identical to the original input data
#    (stored in file "orig.bin"). Check performed for 
#    'test' files only.
#  - The kgds array values from the 'control' and
#    'test' are stored to log files.  The log
#    files must be identical.
#
# If anything fails, the log and binary data are stored
# in a $WORK_DIR subdirectory with a 'failed' extension.
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

    cmp $WORK_CTL/ipsector.namer.bin $WORK_TEST/ipsector.namer.bin
    status=$?
    if ((status != 0))
    then
      echo "N AMERICA IPSECTOR BINARY FILES NOT BIT IDENTICAL. TEST FAILED."
      save_ctl=1
      save_test=1
    fi

    cmp $WORK_CTL/ipsector.no.sect.bin $WORK_TEST/ipsector.no.sect.bin
    status=$?
    if ((status != 0))
    then
      echo "NON-OVERLAP IPSECTOR BINARY FILES NOT BIT IDENTICAL. TEST FAILED."
      save_ctl=1
      save_test=1
    fi

    cmp $WORK_CTL/ipsector.ovlp.sect.bin $WORK_TEST/ipsector.ovlp.sect.bin
    status=$?
    if ((status != 0))
    then
      echo "OVERLAP IPSECTOR BINARY FILES NOT BIT IDENTICAL. TEST FAILED."
      save_ctl=1
      save_test=1
    fi

    cmp $WORK_TEST/orig.bin $WORK_TEST/ipspaste.namer.bin
    status=$?
    if ((status != 0))
    then
      echo "N AMERICA IPSPASTE BINARY FILE AND INPUT FILE NOT BIT IDENTICAL. TEST FAILED."
      save_test=1
    fi

    cmp $WORK_TEST/orig.bin $WORK_TEST/ipspaste.no.sect.bin
    status=$?
    if ((status != 0))
    then
      echo "NON-OVERLAP IPSPASTE BINARY FILE AND INPUT FILE NOT BIT IDENTICAL. TEST FAILED."
      save_test=1
    fi

    cmp $WORK_TEST/orig.bin $WORK_TEST/ipspaste.ovlp.sect.bin
    status=$?
    if ((status != 0))
    then
      echo "OVERLAP IPSPASTE BINARY FILE AND INPUT FILE NOT BIT IDENTICAL. TEST FAILED."
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
      FAILED_DIR=$WORK_CTL/failed.${bytesize}byte.scan${scan}
      mkdir -p $FAILED_DIR
      cp $WORK_CTL/$CTL_LOG  $FAILED_DIR
      cp $WORK_CTL/*.bin     $FAILED_DIR
    fi

    if ((save_test == 1)); then
      FAILED_DIR=$WORK_TEST/failed.${bytesize}byte.scan${scan}
      mkdir -p $FAILED_DIR
      cp $WORK_TEST/$TEST_LOG  $FAILED_DIR
      cp $WORK_TEST/*.bin      $FAILED_DIR
    fi

    rm -f $WORK_CTL/*.bin $WORK_TEST/*.bin

  done  # library byte size loop

  rm -f $WORK_CTL/fort.9 $WORK_TEST/fort.9
  echo

done # scan mode loop

echo
echo "IPSECTOR/IPSPASTE REGRESSION TEST COMPLETED."
echo

exit 0
