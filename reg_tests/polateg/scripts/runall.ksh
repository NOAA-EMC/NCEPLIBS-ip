#!/bin/ksh

#---------------------------------------------------------------
# Test the polateg suite of iplib routines.
#
# Interpolate a global field of temperature to a high-res
# global grid and output the vector gradients to a direct
# access file.  polateg0 - bilinear; polateg1 - bicubic;
# polateg4 - spectral.
#
# Check control and test files for bit-identicalness.
# If not identical, the test is considered "failed".
# The log and binary files are then saved with a ".failed"
# extension.
#---------------------------------------------------------------

#set -x

WORK_DIR=${WORK_DIR:-/stmp/$LOGNAME/regression}

REG_DIR=${REG_DIR:-../..}

EXEC_DIR=$REG_DIR/polateg/exec

INPUT_DATA=$REG_DIR/polateg/data/global_tg3clim.1x1.grb

WORK=$WORK_DIR/polateg
rm -fr $WORK
mkdir -p $WORK
WORK_CTL=${WORK}/ctl
mkdir -p $WORK_CTL
cp $EXEC_DIR/ctl/*.exe $WORK_CTL
cp $INPUT_DATA  $WORK_CTL/fort.9
WORK_TEST=${WORK}/test
mkdir -p $WORK_TEST
cp $EXEC_DIR/test/*.exe $WORK_TEST
cp $INPUT_DATA  $WORK_TEST/fort.9

for bytesize in "4" "8" "d"
do

  save_logs=0

  echo "TEST ${bytesize}-BYTE FLOAT VERSION OF POLATEG ROUTINES"

  cd $WORK_CTL
  polateg_ctl_${bytesize}.exe  > ctl.log
  status=$? 
  if ((status != 0)); then
    echo "CTL RUN FAILED.  CHECK LOG FILE."
    save_logs=1
    if [[ -s core ]];then
       mv core core.${bytesize}.byte
    fi
  fi

  cd $WORK_TEST
  polateg_test_${bytesize}.exe > test.log
  status=$? 
  if ((status != 0)); then
    echo "TEST RUN FAILED.  CHECK LOG FILE"
    save_logs=1
    if [[ -s core ]];then
       mv core core.${bytesize}.byte
    fi
  fi

  cmp $WORK_CTL/polateg0.bin $WORK_TEST/polateg0.bin
  status=$?
  if ((status != 0))
  then
    echo "POLATEG0 BINARY FILES ARE NOT BIT IDENTICAL. TEST FAILED."
    save_logs=1
    if [[ -s $WORK_CTL/polateg0.bin ]];then
      mv $WORK_CTL/polateg0.bin $WORK_CTL/polateg0.${bytesize}.byte.bin.failed
    else
      echo "CTL POLATEG0 BINARY FILE MISSING."
    fi
    if [[ -s $WORK_TEST/polateg0.bin ]];then
      mv $WORK_TEST/polateg0.bin $WORK_TEST/polateg0.${bytesize}.byte.bin.failed
    else
      echo "TEST POLATEG0 BINARY FILE MISSING."
    fi
  else
    rm -f $WORK_CTL/polateg0.bin $WORK_TEST/polateg0.bin
  fi

  cmp $WORK_CTL/polateg1.bin $WORK_TEST/polateg1.bin
  status=$?
  if ((status != 0))
  then
    echo "POLATEG1 BINARY FILES ARE NOT BIT IDENTICAL. TEST FAILED."
    save_logs=1
    if [[ -s $WORK_CTL/polateg1.bin ]];then
      mv $WORK_CTL/polateg1.bin $WORK_CTL/polateg1.${bytesize}.byte.bin.failed
    else
      echo "CTL POLATEG1 BINARY FILE MISSING."
    fi
    if [[ -s $WORK_TEST/polateg1.bin ]];then
      mv $WORK_TEST/polateg1.bin $WORK_TEST/polateg1.${bytesize}.byte.bin.failed
    else
      echo "TEST POLATEG1 BINARY FILE MISSING."
    fi
  else
    rm -f $WORK_CTL/polateg1.bin $WORK_TEST/polateg1.bin
  fi

  cmp $WORK_CTL/polateg4.bin $WORK_TEST/polateg4.bin
  status=$?
  if ((status != 0))
  then
    echo "POLATEG4 BINARY FILES ARE NOT BIT IDENTICAL. TEST FAILED."
    save_logs=1
    if [[ -s $WORK_CTL/polateg4.bin ]];then
      mv $WORK_CTL/polateg4.bin $WORK_CTL/polateg4.${bytesize}.byte.bin.failed
    else
      echo "CTL POLATEG4 BINARY FILE MISSING."
    fi
    if [[ -s $WORK_TEST/polateg4.bin ]];then
      mv $WORK_TEST/polateg4.bin $WORK_TEST/polateg4.${bytesize}.byte.bin.failed
    else
      echo "TEST POLATEG4 BINARY FILE MISSING."
    fi
  else
    rm -f $WORK_CTL/polateg4.bin $WORK_TEST/polateg4.bin
  fi

  if ((save_logs == 1));then
    mv $WORK_CTL/ctl.log  $WORK_CTL/ctl.${bytesize}.byte.log.failed
    mv $WORK_TEST/test.log  $WORK_TEST/test.${bytesize}.byte.log.failed
  else
    rm -f $WORK_CTL/ctl.log $WORK_TEST/test.log
  fi

done

echo
echo "TESTS COMPLETED."

exit 0
