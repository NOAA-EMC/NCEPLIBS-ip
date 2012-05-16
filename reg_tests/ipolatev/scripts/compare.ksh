#!/bin/ksh

# ensure single and multiple thread files are bit identical

#set -x

echo
echo "ENSURE SINGLE AND MULTIPLE THREAD IPOLATEV FILES ARE BIT IDENTICAL"
echo

typeset -L1 machine
machine=$(hostname)

if [[ $machine == "f" ]];then    # zeus
  export WORK_DIR=${WORK_DIR:-/scratch2/portfolios/NCEPDEV/stmp/$LOGNAME/regression}
else   # cirrus/stratus
  export WORK_DIR=${WORK_DIR:-/stmp/$LOGNAME/regression}
fi

WORK1=$WORK_DIR/ipolatev.1threads/test
WORK2=$WORK_DIR/ipolatev.4threads/test

cd $WORK1

for binfile in *.bin
do
  echo "COMPARE FILE " $binfile
  cmp $binfile  $WORK2/$binfile
  status=$?
  if ((status != 0))
  then
    echo $binfile "NOT BIT IDENTIAL. TEST FAILED."
  fi
done

echo
echo TEST COMPLETED
echo

exit 0
