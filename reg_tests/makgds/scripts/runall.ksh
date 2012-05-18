#!/bin/ksh

#set -x

echo
echo BEGIN MAKGDS REGRESSION TEST
echo
 
typeset -L4 machine
machine=$(dnsdomainname)

if [[ $machine == "zeus" ]];then    # zeus
  WORK_DIR=${WORK_DIR:-/scratch2/portfolios/NCEPDEV/stmp/$LOGNAME/regression}
else   # cirrus/stratus
  WORK_DIR=${WORK_DIR:-/stmp/$LOGNAME/regression}
fi

REG_DIR=${REG_DIR:-../..}

EXEC_DIR=$REG_DIR/makgds/exec

WORK=$WORK_DIR/makgds
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
  save_ctl_log=0
  save_test_log=0
  echo TEST ${bytesize}-BYTE FLOAT VERSION OF ROUTINE MAKGDS
  cd $WORK_CTL
  makgds_ctl_${bytesize}.exe  > ctl.log
  grep -Eq 'BAD|ERROR' $WORK_CTL/ctl.log
  status=$?
  if ((status == 0)); then
    echo PROBLEM WITH CTL RUN. CHECK LOG FILE.
    save_ctl_log=1
  fi
  cd $WORK_TEST
  makgds_test_${bytesize}.exe > test.log
  grep -Eq 'BAD|ERROR' $WORK_TEST/test.log
  status=$?
  if ((status == 0)); then
    echo PROBLEM WITH TEST RUN. CHECK LOG FILE.
    save_test_log=1
  fi
  cmp $WORK_CTL/ctl.log $WORK_TEST/test.log
  status=$?
  if ((status != 0))
  then
    echo LOG FILES NOT BIT IDENTIAL. TEST FAILED.
    save_ctl_log=1
    save_test_log=1
  fi

  if ((save_ctl_log == 1)); then
    mv $WORK_CTL/ctl.log $WORK_CTL/ctl.${bytesize}byte.log.failed
  fi

  if ((save_test_log == 1)); then
    mv $WORK_TEST/test.log $WORK_TEST/test.${bytesize}byte.log.failed
  fi

done

echo
echo MAKGDS REGRESSION TEST COMPLETED.
echo

exit 0
