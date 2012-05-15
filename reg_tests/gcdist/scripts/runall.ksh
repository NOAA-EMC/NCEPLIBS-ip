#!/bin/ksh

#set -x

echo
echo BEGIN GCDIST/MOVECT REGRESSION TEST
echo

typeset -L1 machine
machine=$(hostname)

if [[ $machine == "f" ]];then    # zeus
  export REG_DIR=${REG_DIR:-../..}
  export WORK_DIR=${WORK_DIR:-/scratch2/portfolios/NCEPDEV/stmp/$LOGNAME}
else   # cirrus/stratus
  export REG_DIR=${REG_DIR:-../..}
  export WORK_DIR=${WORK_DIR:-/stmp/$LOGNAME/regression}
fi

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

for bytesize in "4" "8" "d"
do
  echo TEST ${bytesize}-BYTE FLOAT VERSION OF ROUTINES GCDIST/MOVECT
  cd $WORK_CTL
  CTL_LOG=ctl.${bytesize}byte.log
  gcdist_ctl_${bytesize}.exe  > $CTL_LOG
  cd $WORK_TEST
  TEST_LOG=test.${bytesize}byte.log
  gcdist_test_${bytesize}.exe > $TEST_LOG
  cmp $WORK_CTL/$CTL_LOG $WORK_TEST/$TEST_LOG
  status=$?
  if ((status != 0))
  then
    echo LOG FILES NOT BIT IDENTIAL. TEST FAILED.
    echo CHECK LOG FILES SAVED IN WORK DIRECTORY.
    mv $WORK_CTL/$CTL_LOG $WORK_CTL/${CTL_LOG}.failed
    mv $WORK_TEST/$TEST_LOG $WORK_TEST/${TEST_LOG}.failed
  fi
done

echo
echo GCDIST/MOVECT REGRESSION TEST COMPLETED
echo

exit 0
