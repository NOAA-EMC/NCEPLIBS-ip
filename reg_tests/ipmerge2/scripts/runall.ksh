#!/bin/ksh

#set -x

echo
echo BEGIN IPMERGE2 REGRESSION TEST
echo

typeset -L4 machine
machine=$(dnsdomainname)

REG_DIR=${REG_DIR:-../..}

if [[ $machine == "zeus" ]];then    # zeus
  WORK_DIR=${WORK_DIR:-/scratch2/portfolios/NCEPDEV/stmp/$LOGNAME/regression}
else   # cirrus/stratus
  WORK_DIR=${WORK_DIR:-/stmp/$LOGNAME/regression}
fi

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

for bytesize in "4" "8" "d"
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
  fi
done

echo
echo IPMERGE2 REGRESSION TEST COMPLETED.
echo

exit 0
