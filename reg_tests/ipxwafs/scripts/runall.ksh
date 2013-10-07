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
# then the regression test has failed.  when this happens, the output file is
# saved in a $WORK_DIR subdirectory with a "failed" extension.
#-------------------------------------------------------------------------------

#set -x

echo
echo BEGIN IPXWAFS REGRESSION TEST
echo

WORK_DIR=${WORK_DIR:-/stmp/$LOGNAME/regression}

REG_DIR=${REG_DIR:-../..}

TEST_EXEC_DIR=$REG_DIR/ipxwafs/exec
CTL_EXEC_DIR=$REG_DIR/ipxwafs/exec

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

reg_test_failed=0

echo
echo CONVERT FROM REGULAR GRID TO WAFS GRIDS.
for bytesize in "4" "8" "d"  # test all byte versions of iplib
do
  echo
  echo TEST $bytesize BYTE VERSION OF LIBRARY.
  echo
  for grid in "37" "38" "39" "40" "41" "42" "43" "44"
  do
    echo CONVERT TO WAFS GRID ${grid}.
    for ipopt in "0" 
    do
      echo TEST INTERPOLATION OPTION $ipopt

      ctl_failed=0
      test_failed=0

      cd $WORK_TEST
      copygb_test_${bytesize} -g${grid} -i${ipopt} -x ../input_data/600mb.temp.grb ${grid}.grb 
      status=$?
      if ((status != 0))
      then
        echo "** TEST RUN FAILED **"
        reg_test_failed=1
        test_failed=1
      fi

      cd $WORK_CTL
      copygb_ctl_${bytesize} -g${grid} -i${ipopt} -x ../input_data/600mb.temp.grb ${grid}.grb
      status=$?
      if ((status != 0))
      then
        echo "** CONTROL RUN FAILED **"
        reg_test_failed=1
        ctl_failed=1
      fi

      if ((ctl_failed == 0 && test_failed == 0));then
        cmp $WORK_CTL/${grid}.grb $WORK_TEST/${grid}.grb
        status=$?
        if ((status != 0))
        then
          echo "** GRIB FILES NOT BIT IDENTICAL. REGRESSION TEST FAILED."
          FAILED_DIR=$WORK_CTL/failed.regl.to.wafs
          mkdir -p $FAILED_DIR
          mv $WORK_CTL/${grid}.grb $FAILED_DIR/${grid}.${bytesize}byte.grb
          FAILED_DIR=$WORK_TEST/failed.regl.to.wafs
          mkdir -p $FAILED_DIR
          mv $WORK_TEST/${grid}.grb $FAILED_DIR/${grid}.${bytesize}byte.grb
          reg_test_failed=1
        fi
      fi

      rm -f $WORK_CTL/${grid}.grb $WORK_TEST/${grid}.grb

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

      ctl_failed=0
      test_failed=0

      cd $WORK_TEST
      copygb_test_${bytesize} -g3 -i${ipopt} -x ../input_data/wafs.${grid}.grb  grid3.from${grid}.grb
      status=$?
      if ((status != 0))
      then
        echo "** TEST RUN FAILED **"
        reg_test_failed=1
        test_failed=1
      fi

      cd $WORK_CTL
      copygb_ctl_${bytesize} -g3 -i${ipopt} -x ../input_data/wafs.${grid}.grb  grid3.from${grid}.grb
      status=$?
      if ((status != 0))
      then
        echo "** CONTROL RUN FAILED **"
        reg_test_failed=1
        ctl_failed=1
      fi

      if ((ctl_failed == 0 && test_failed == 0));then
        cmp $WORK_CTL/grid3.from${grid}.grb $WORK_TEST/grid3.from${grid}.grb
        status=$?
        if ((status != 0))
        then
          echo "** GRIB FILES NOT BIT IDENTICAL. REGRESSION TEST FAILED."
          FAILED_DIR=$WORK_CTL/failed.wafs.to.regl
          mkdir -p $FAILED_DIR
          mv $WORK_CTL/grid3.from${grid}.grb $FAILED_DIR/grid3.from${grid}.${bytesize}byte.grb
          FAILED_DIR=$WORK_TEST/failed.wafs.to.regl
          mkdir -p $FAILED_DIR
          mv $WORK_TEST/grid3.from${grid}.grb $FAILED_DIR/grid3.from${grid}.${bytesize}byte.grb
          reg_test_failed=1
        fi
      fi

      rm -f $WORK_CTL/grid3.from${grid}.grb $WORK_TEST/grid3.from${grid}.grb

    done
  done
done

if ((reg_test_failed == 0));then
  echo
  echo "<<< IPXWAFS REGRESSION TEST PASSED. >>>"
  echo
else
  echo
  echo "<<< IPXWAFS REGRESSION TEST FAILED. >>>"
  echo
fi

exit 0
