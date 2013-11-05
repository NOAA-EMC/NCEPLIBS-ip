#!/bin/ksh

#--------------------------------------------------------------------
# Part 2 of the ipolatev regression test.
#
# Ensure single and multiple thread files are bit identical.
# If not, the regression test fails.
#--------------------------------------------------------------------

#set -x

echo
echo "ENSURE SINGLE AND MULTIPLE THREAD IPOLATEV FILES ARE BIT IDENTICAL"
echo

WORK_DIR=${WORK_DIR:-/stmp/$LOGNAME/regression}

WORK1=$WORK_DIR/ipolatev.1threads/test
WORK2=$WORK_DIR/ipolatev.4threads/test

cd $WORK1

failed=0

for binfile in *.bin
do
  echo "COMPARE FILE " $binfile
  cmp $binfile  $WORK2/$binfile
  status=$?
  if ((status != 0))
  then
    echo $binfile "NOT BIT IDENTICAL. TEST FAILED."
    failed=1
  fi
done

if ((failed == 0));then
  echo
  echo "<<< IPOLATEV FILE THREAD TEST PASSED. >>>"
  echo
else
  echo
  echo "<<< IPOLATEV FILE THREAD TEST FAILED. >>>"
  echo
fi

exit 0
