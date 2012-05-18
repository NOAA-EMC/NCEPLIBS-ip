#!/bin/ksh

#set -x

if (($# > 0))
then
  num_threads=$1
  echo
  echo "BEGIN IPOLATES REGRESSION TEST WITH " $num_threads "THREADS"
else
  echo "ENTER NUMBER OF THREADS"
  exit 99
fi

typeset -L4 machine
machine=$(dnsdomainname)

if [[ $machine == "zeus" ]];then    # zeus
  WORK_DIR=${WORK_DIR:-/scratch2/portfolios/NCEPDEV/stmp/$LOGNAME/regression}
else   # cirrus/stratus
  WORK_DIR=${WORK_DIR:-/stmp/$LOGNAME/regression}
fi

REG_DIR=${REG_DIR:-../..}

EXEC_DIR=$REG_DIR/ipolates/exec
INPUT_DATA=$REG_DIR/ipolates/data/global_tg3clim.1x1.grb

WORK=$WORK_DIR/ipolates.${num_threads}threads
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

for grids in "3" "8" "127" "203" "205" "212" "218" 
do
  echo
  for option in "0" "1" "2" "3" "4" "6"  # interpolation option
  do
    for bytesize in "4" "8" "d"
    do
      echo TEST ${bytesize}-BYTE VERSION FOR GRID $grids AND INTERP OPTION $option
      cd $WORK_CTL
      ipolates_ctl_${bytesize}.exe "$grids" "$option" > ctl.log
      cd $WORK_TEST
      ipolates_test_${bytesize}.exe "$grids" "$option" > test.log

      save_ctl_log=0
      save_test_log=0

      cmp $WORK_CTL/grid${grids}.opt${option}.bin $WORK_TEST/grid${grids}.opt${option}.bin
      status=$?
      if ((status != 0))
      then
        echo BINARY FILES NOT BIT IDENTIAL. TEST FAILED.
        mv $WORK_CTL/grid${grids}.opt${option}.bin $WORK_CTL/grid${grids}.opt${option}.${bytesize}byte.bin.failed
        mv $WORK_TEST/grid${grids}.opt${option}.bin $WORK_TEST/grid${grids}.opt${option}.${bytesize}byte.bin.failed
        save_ctl_log=1
        save_test_log=1
      else
        mv $WORK_CTL/grid${grids}.opt${option}.bin $WORK_CTL/grid${grids}.opt${option}.${bytesize}byte.bin
        mv $WORK_TEST/grid${grids}.opt${option}.bin $WORK_TEST/grid${grids}.opt${option}.${bytesize}byte.bin
      fi

      grep -Eq 'BAD|ERROR' $WORK_CTL/ctl.log
      status=$?
      if ((status == 0));then
        echo PROBLEM WITH CTL RUN. CHECK LOG FILE.
        save_ctl_log=1
      fi

      grep -Eq 'BAD|ERROR' $WORK_TEST/test.log
      status=$?
      if ((status == 0));then
        echo PROBLEM WITH TEST RUN. CHECK LOG FILE.
        save_test_log=1
      fi

      if ((save_ctl_log == 1));then
        mv $WORK_CTL/ctl.log $WORK_CTL/ctl.grid${grids}.opt${option}.${bytesize}byte.log.failed
      else
        rm -f $WORK_CTL/ctl.log
      fi

      if ((save_test_log == 1));then
        mv $WORK_TEST/test.log $WORK_TEST/test.grid${grids}.opt${option}.${bytesize}byte.log.failed
      else
        rm -f $WORK_TEST/test.log
      fi

    done
  done
done

echo
echo "IPOLATES REGRESSION TEST WITH " $num_threads "THREADS COMPLETED."
echo

exit 0
