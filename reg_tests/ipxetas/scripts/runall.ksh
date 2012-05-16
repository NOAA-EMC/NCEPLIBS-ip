#!/bin/ksh

#------------------------------------------------------------------------
# test iplib routine ipxetas.  
#
# the routine is called four times:
#
# 1) create a staggered mass grid from the full grid.
# 2) create a staggered velocity grid from the full grid.
# 3) create a full grid from the staggered mass grid created by step (1)
# 4) create a full grid from the staggered vel grid created by step (2)
#
# output from steps (1) and (2) is in the file named "staggered.bin"
# output from steps (3) and (4) is in the file named "full.bin"
#
#------------------------------------------------------------------------

#set -x

echo
echo BEGIN REGRESSION TEST FOR IPXETAS
echo

typeset -L1 machine
machine=$(hostname)

if [[ $machine == "f" ]];then    # zeus
  export REG_DIR=${REG_DIR:-../..}
  export WORK_DIR=${WORK_DIR:-/scratch2/portfolios/NCEPDEV/stmp/$LOGNAME/regression}
else   # cirrus/stratus
  export REG_DIR=${REG_DIR:-../..}
  export WORK_DIR=${WORK_DIR:-/stmp/$LOGNAME/regression}
fi

EXEC_DIR=${REG_DIR}/ipxetas/exec
INPUT_DATA=$REG_DIR/ipxetas/data/green.202.grb

WORK=${WORK_DIR}/ipxetas
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

for bytesize in "4" "8" "d"
do

  echo TEST ${bytesize}-BYTE VERSION OF IPXETAS
  cd $WORK_CTL
  CTL_LOG=ctl.${bytesize}byte.log
  ipxetas_ctl_${bytesize}.exe  > $CTL_LOG
  cd $WORK_TEST
  TEST_LOG=test.${bytesize}byte.log
  ipxetas_test_${bytesize}.exe > $TEST_LOG

  save_ctl_log=0
  save_test_log=0

  cmp $WORK_CTL/staggered.bin $WORK_TEST/staggered.bin
  status=$?
  if ((status != 0))
  then
    echo STAGGERED GRID BINARY FILES NOT BIT IDENTICAL. TEST FAILED.
    mv $WORK_CTL/staggered.bin $WORK_CTL/staggered.${bytesize}byte.bin.failed
    mv $WORK_TEST/staggered.bin $WORK_TEST/staggered.${bytesize}byte.bin.failed
    save_ctl_log=1
    save_test_log=1
  fi

  cmp $WORK_CTL/full.bin $WORK_TEST/full.bin
  status=$?
  if ((status != 0))
  then
    echo FULL GRID BINARY FILES NOT BIT IDENTICAL. TEST FAILED.
    mv $WORK_CTL/full.bin $WORK_CTL/full.${bytesize}byte.bin.failed
    mv $WORK_TEST/full.bin $WORK_TEST/full.${bytesize}byte.bin.failed
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
    mv $WORK_CTL/$CTL_LOG  $WORK_CTL/${CTL_LOG}.failed
  fi

  if ((save_test_log == 1)); then
    mv $WORK_TEST/$TEST_LOG  $WORK_TEST/${TEST_LOG}.failed
  fi

  rm -f $WORK_CTL/staggered.bin $WORK_TEST/staggered.bin
  rm -f $WORK_CTL/full.bin $WORK_TEST/full.bin

done

echo
echo IPXETAS REGRESSION TEST COMPLETED.
echo

exit 0
