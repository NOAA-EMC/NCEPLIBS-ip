#!/bin/ksh

#set -x

echo
echo BEGIN GAUSSLAT REGRESSION TEST
echo

typeset -L4 machine
machine=$(dnsdomainname)

REG_DIR=${REG_DIR:-../..}

if [[ $machine == "zeus" ]];then    # zeus
  WORK_DIR=${WORK_DIR:-/scratch2/portfolios/NCEPDEV/stmp/$LOGNAME/regression}
else   # cirrus/stratus
  WORK_DIR=${WORK_DIR:-/stmp/$LOGNAME/regression}
fi 

EXEC_DIR=$REG_DIR/gausslat/exec

WORK=${WORK_DIR}/gausslat
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
  echo TEST ${bytesize}-BYTE FLOAT VERSION OF ROUTINE GAUSSLAT
  cd $WORK_CTL
  CTL_LOG=ctl.${bytesize}byte.log
  gausslat_ctl_${bytesize}.exe  > $CTL_LOG
  cd $WORK_TEST
  TEST_LOG=test.${bytesize}byte.log
  gausslat_test_${bytesize}.exe > $TEST_LOG
  cmp $WORK_CTL/$CTL_LOG $WORK_TEST/$TEST_LOG
  status=$?
  if ((status != 0))
  then
    echo LOG FILES NOT BIT IDENTIAL. TEST FAILED.
    echo CHECK LOG FILE SAVED IN WORK DIRECTORY.
    mv $WORK_CTL/$CTL_LOG $WORK_CTL/${CTL_LOG}.failed
    mv $WORK_TEST/$TEST_LOG $WORK_TEST/${TEST_LOG}.failed
  fi
done

echo
echo GAUSSLAT REGRESSION TEST COMPLETED.
echo

exit 0
