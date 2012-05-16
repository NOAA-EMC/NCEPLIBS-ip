#!/bin/ksh

# test the gdswiz and gdswzd suite of routines.

#set -x

echo
echo BEGIN GDSWIZ/WZD REGRESSION TEST
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

EXEC_DIR=$REG_DIR/gdswiz_wzd/exec

WORK=${WORK_DIR}/gdswiz_wzd
rm -fr $WORK
mkdir -p $WORK
WORK_CTL=${WORK}/ctl
mkdir -p $WORK_CTL
cp $EXEC_DIR/ctl/*exe $WORK_CTL
WORK_TEST=${WORK}/test
mkdir -p $WORK_TEST
cp $EXEC_DIR/test/*exe $WORK_TEST

for routine in "WIZ" "WZD"
do
  echo
  echo RUN REGRESSION TEST FOR GDS${routine} ROUTINES 
  for grids in "3" "8" "203" "127" "212" "213" "218" "205" "201" "202" "222"
  do
    echo
    for bytesize in "4" "8" "d"
    do

      echo TEST GRID $grids ${bytesize}-BYTE FLOAT VERSION
      cd $WORK_CTL
      CTL_LOG=ctl.${routine}.${bytesize}byte.grid${grids}.log
      gdswiz_wzd_ctl_${bytesize}.exe "$routine" "$grids" > $CTL_LOG
      cd $WORK_TEST
      TEST_LOG=test.${routine}.${bytesize}byte.grid${grids}.log
      gdswiz_wzd_test_${bytesize}.exe "$routine" "$grids" > $TEST_LOG

      save_ctl_log=0
      save_test_log=0

      cmp $WORK_CTL/grid${grids}.bin $WORK_TEST/grid${grids}.bin
      status=$?
      if ((status != 0))
      then
        echo "** BINARY FILES NOT BIT IDENTICAL. TEST FAILED."
        mv $WORK_CTL/grid${grids}.bin $WORK_CTL/grid${grids}.${routine}.${bytesize}byte.bin.failed
        mv $WORK_TEST/grid${grids}.bin $WORK_TEST/grid${grids}.${routine}.${bytesize}byte.bin.failed
        save_ctl_log=1
        save_test_log=1
      fi

      grep -Eq 'BAD|ERROR' $WORK_CTL/$CTL_LOG
      status=$?
      if ((status == 0)); then
        echo "** PROBLEM WITH CTL RUN. CHECK LOG FILE."
        save_ctl_log=1
      fi

      grep -Eq 'BAD|ERROR' $WORK_TEST/$TEST_LOG
      status=$?
      if ((status == 0)); then
        echo "** PROBLEM WITH TEST RUN. CHECK LOG FILE."
        save_test_log=1
      fi

      if ((save_ctl_log == 1));then
        mv $WORK_CTL/$CTL_LOG $WORK_CTL/${CTL_LOG}.failed
      fi

      if ((save_test_log == 1));then
         mv $WORK_TEST/$TEST_LOG $WORK_TEST/${TEST_LOG}.failed
      fi

      rm -f $WORK_CTL/grid${grids}.bin $WORK_TEST/grid${grids}.bin
      rm -f $WORK_CTL/$CTL_LOG $WORK_TEST/$TEST_LOG

    done
  done
done

echo
echo GDSWIZ/WZD REGRESSION TEST COMPLETED.
echo

exit 0
