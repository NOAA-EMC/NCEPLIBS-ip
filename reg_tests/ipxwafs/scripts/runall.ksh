#!/bin/ksh

#-------------------------------------------------------------------------------
# test ip routine ipxwafs by transforming a global grid of
# scalars (on ncep grid 3) to wafs grids 37 thru 44 using a specially 
# modified version of copygb.
#
# after the global to wafs grid transforms are completed, copygb is 
# invoked again to go from each wafs grid back to ncep grid 3.
#
# note: routine ipxwafs and ipxwafs2 are the same except the latter
# accounts for bitmaps.  the ops version of copygb uses the latter
# routine for all interpolation options except '2' (neighbor).
# the specially modified version calls ipxwafs instead.  since the
# test dataset does not have a bitmap (a field of 600 mb temperatures)
# the specially modified copygb and ops copygb give the same answer.
#
# if the output files from the control and test are not bit identical,
# then the test has failed.  if a test fails, the output file is saved
# in WORK_DIR with a ".failed" extension.
#-------------------------------------------------------------------------------

#set -x

echo
echo BEGIN IPXWAFS REGRESSION TEST
echo

WORK_DIR=${WORK_DIR:-/stmp/$LOGNAME/regression}

REG_DIR=${REG_DIR:-../..}

TEST_EXEC_DIR=$REG_DIR/ipxwafs/exec/test
CTL_EXEC_DIR=$REG_DIR/ipxwafs/exec/ctl

DATA_DIR=$REG_DIR/ipxwafs2_3/data

WORK=$WORK_DIR/ipxwafs
rm -fr $WORK
mkdir -p $WORK
mkdir -p $WORK/input_data
cp $DATA_DIR/* $WORK/input_data
WORK_CTL=${WORK}/ctl
mkdir -p $WORK_CTL
cp $CTL_EXEC_DIR/copygb_ctl* $WORK_CTL
WORK_TEST=${WORK}/test
mkdir -p $WORK_TEST
cp $TEST_EXEC_DIR/copygb_test* $WORK_TEST

echo
echo CONVERT FROM REGULAR GRID TO WAFS GRIDS.
for bytesize in "4" "8" "d"  # test all byte versions of iplib
do
  echo
  echo TEST $bytesize VERSION OF LIBRARY.
  echo
  for grid in "37" "38" "39" "40" "41" "42" "43" "44"
  do
    echo CONVERT TO WAFS GRID ${grid}.
    for ipopt in "0" 
    do
      echo TEST INTERPOLATION OPTION $ipopt

      cd $WORK_TEST
      copygb_test_${bytesize} -g${grid} -i${ipopt} -x ../input_data/600mb.temp.grb ${grid}.grb 
      status=$?
      if ((status != 0))
      then
        echo "** TEST FAILED **"
        continue
      fi

      cd $WORK_CTL
      copygb_ctl_${bytesize} -g${grid} -i${ipopt} -x ../input_data/600mb.temp.grb ${grid}.grb
      status=$?
      if ((status != 0))
      then
        echo "** TEST FAILED **"
        continue
      fi

      cmp $WORK_CTL/${grid}.grb $WORK_TEST/${grid}.grb
      status=$?
      if ((status != 0))
      then
        echo "** GRIB FILES NOT BIT IDENTICAL. TEST FAILED."
        mv $WORK_CTL/${grid}.grb $WORK_CTL/${grid}.${bytesize}byte.grb.failed
        mv $WORK_TEST/${grid}.grb $WORK_TEST/${grid}.${bytesize}byte.grb.failed
      else
        rm -f $WORK_CTL/${grid}.grb $WORK_TEST/${grid}.grb
      fi

    done
  done
done

echo
echo CONVERT FROM WAFS GRIDS TO REGULAR GRIDS.
for bytesize in "4" "8" "d"
do
  echo
  echo TEST $bytesize BYTE VERSION OF LIBRARY.
  echo
  for grid in "37" "38" "39" "40" "41" "42" "43" "44"
  do
    echo CONVERT FROM WAFS GRID ${grid}.
    for ipopt in "0" 
    do
      echo TEST INTERPOLATION OPTION $ipopt

      cd $WORK_TEST
      copygb_test_${bytesize} -g3 -i${ipopt} -x ../input_data/wafs.${grid}.grb  grid3.grb
      status=$?
      if ((status != 0))
      then
        echo "** TEST FAILED **"
        continue
      fi

      cd $WORK_CTL
      copygb_ctl_${bytesize} -g3 -i${ipopt} -x ../input_data/wafs.${grid}.grb  grid3.grb
      status=$?
      if ((status != 0))
      then
        echo "** TEST FAILED **"
        continue
      fi

      cmp $WORK_CTL/grid3.grb $WORK_TEST/grid3.grb
      status=$?
      if ((status != 0))
      then
        echo "** GRIB FILES NOT BIT IDENTICAL. TEST FAILED."
      fi
      rm -f $WORK_CTL/grid3.grb $WORK_TEST/grid3.grb

    done
  done
done

echo
echo IPXWAFS REGRESSION TEST COMPLETED.
echo

exit 0
