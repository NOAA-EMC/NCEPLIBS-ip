#!/bin/ksh

#-------------------------------------------------------------------------------
# Test ip2 routines ipxwafs, ipxwafs2 and ipxwafs3 by transforming between
# full and thinned WAFS grids.  The full grids have the same number of
# points in each row.  In the thinned grids, the number of points in each
# row decrease toward the pole. These transforms are performed using a
# Fortran program.
#
# The program reads in a grib 2 file of data on the full/thinned WAFS
# grid and transforms it to its thinned/full counterpart.  Routine
# ipxwafs transforms data without bitmaps.  Routines ipxwafs2 and
# ipxwafs3 transform data with bitmaps - the former uses bilinear
# interpolation and the latter uses nearest neighbor.  The transformed
# data is written to a grib 2 file.
#
# The program executables are located under the ./exec subdirectory.
# There is one executable for all three byte versions of the
# 'control' and 'test' ip2 library:
#
# The three byte versions of the library are:
#  > 4 byte integer/4 byte float  ($bytesize=4)
#  > 8 byte integer/8 byte float  ($bytesize=8)
#  > 8 byte float/4 byte integer  ($bytesize=d)
#
# The input data for the program is located in the /data subdirectory.
# All input data are in grib 2 format.  The files are:
#
#   wafs.37.full.bitmap.grb2 (soil T on full WAFS grid #37 - land only)
#   wafs.37.full.grb2        (600 mb T on full WAFS grid #37)
#   wafs.37.thin.bitmap.grb2 (soil T on thinned WAFS grid #37 - land only)
#   wafs.37.thin.grb2        (600 mb T on thinned WAFS grid #37)
#   wafs.44.full.bitmap.grb2 (soil T on full WAFS grid #44 - land only)
#   wafs.44.full.grb2        (600 mb T on full WAFS grid #44)
#   wafs.44.thin.bitmap.grb2 (soil T on thinned WAFS grid #44 - land only)
#   wafs.44.thin.grb2        (600 mb T on thinned WAFS grid #44)
#
# The 'control' and 'test' executables run in their own working directories:
# $WORK_CTL and $WORK_TEST.
#
# If the output files (grib 2 format) from the control and test libraries
# are not bit identical, then the regression test has failed.  When this
# happens, the log and output files are saved in the working directories
# with the following naming convention:
#
#  > log.${transform}.${routine}.grid${wafs_grid}.${bytesize}byte.failed
#  > output.${transform}.${routine}.grid${wafs_grid}.${bytesize}byte.failed.grb2
#
# where: $tranform is "thin2full" or "full2thin"
# where: $routine is "ipxwafs", "ipxwafs2" or "ipxwafs3"
# where: $wafs_grid is "37" or "44"
# where: $bytesize is the library byte size: "4" "8" "d"
#
# This script is run by the "Runall.${machine}" driver script located
# in /reg_tests.
#-------------------------------------------------------------------------------

#set -x

echo
echo BEGIN IPXWAFS REGRESSION TEST
echo

WORK_DIR=${WORK_DIR:-/stmpp1/$LOGNAME/regression}

REG_DIR=${REG_DIR:-../..}
EXEC_DIR=$REG_DIR/ipxwafs/exec
DATA_DIR=$REG_DIR/ipxwafs/data

WORK=$WORK_DIR/ipxwafs
rm -fr $WORK
mkdir -p $WORK

WORK_CTL=${WORK}/ctl
mkdir -p $WORK_CTL
cp $EXEC_DIR/ipxwafs_ctl* $WORK_CTL

WORK_TEST=${WORK}/test
mkdir -p $WORK_TEST
cp $EXEC_DIR/ipxwafs_test* $WORK_TEST

reg_test_failed=0

echo CONVERT FROM THINNED WAFS GRIDS TO FULL GRIDS
for bytesize in "4" "8" "d"  # test all byte versions of library
do

  echo
  echo TEST $bytesize BYTE VERSION OF LIBRARY.

  for wafs_grid in "37" "44"
  do

    echo
    echo CONVERT FROM WAFS GRID ${wafs_grid}.
    echo
    for routine in "ipxwafs" "ipxwafs2" "ipxwafs3"
    do

      echo TEST ROUTINE $routine

      case $routine in
      "ipxwafs")
         INPUT_DATA=$DATA_DIR/wafs.${wafs_grid}.thin.grb2
         cp $INPUT_DATA $WORK_CTL/fort.9
         cp $INPUT_DATA $WORK_TEST/fort.9
         opt=1 ;;
      "ipxwafs2")
         INPUT_DATA=$DATA_DIR/wafs.${wafs_grid}.thin.bitmap.grb2
         cp $INPUT_DATA $WORK_CTL/fort.9
         cp $INPUT_DATA $WORK_TEST/fort.9
         opt=2 ;;
      "ipxwafs3")
         INPUT_DATA=$DATA_DIR/wafs.${wafs_grid}.thin.bitmap.grb2
         cp $INPUT_DATA $WORK_CTL/fort.9
         cp $INPUT_DATA $WORK_TEST/fort.9
         opt=3 ;;
      esac

      ctl_failed=0
      test_failed=0

      cd $WORK_TEST
      ./ipxwafs_test_${bytesize} ${opt} > log
      status=$?
      if ((status != 0))
      then
        echo "** TEST RUN FAILED **"
        mv log log.thin2full.${routine}.grid${wafs_grid}.${bytesize}byte.failed
        reg_test_failed=1
        test_failed=1
      fi

      cd $WORK_CTL
      ./ipxwafs_ctl_${bytesize} ${opt} > log
      status=$?
      if ((status != 0))
      then
        echo "** CONTROL RUN FAILED **"
        mv log log.thin2full.${routine}.grid${wafs_grid}.${bytesize}byte.failed
        reg_test_failed=1
        ctl_failed=1
      fi

      if ((ctl_failed == 0 && test_failed == 0));then
        cmp $WORK_CTL/output.grb2 $WORK_TEST/output.grb2
        status=$?
        if ((status != 0))
        then
          echo "** GRIB FILES NOT BIT IDENTICAL. REGRESSION TEST FAILED."
          mv $WORK_CTL/log $WORK_CTL/log.thin2full.${routine}.grid${wafs_grid}.${bytesize}byte.failed
          mv $WORK_CTL/output.grb2 $WORK_CTL/output.thin2full.${routine}.grid${wafs_grid}.${bytesize}byte.failed.grb2
          mv $WORK_TEST/log $WORK_TEST/log.thin2full.${routine}.grid${wafs_grid}.${bytesize}byte.failed
          mv $WORK_TEST/output.grb2 $WORK_TEST/output.thin2full.${routine}.grid${wafs_grid}.${bytesize}byte.failed.grb2
          reg_test_failed=1
        fi
      fi

      rm -f $WORK_CTL/output.grb2
      rm -f $WORK_CTL/log
      rm -f $WORK_CTL/fort.9
      rm -f $WORK_TEST/output.grb2
      rm -f $WORK_TEST/log
      rm -f $WORK_TEST/fort.9

    done
  done
done

echo
echo CONVERT FROM FULL WAFS GRIDS TO THINNED GRIDS

for bytesize in "4" "8" "d"  # test all byte versions of library
do

  echo
  echo TEST $bytesize BYTE VERSION OF LIBRARY.

  for wafs_grid in "37" "44"
  do

    echo
    echo CONVERT TO WAFS GRID ${wafs_grid}.
    echo
    for routine in "ipxwafs" "ipxwafs2" "ipxwafs3"
    do

      echo TEST ROUTINE $routine

      case $routine in
      "ipxwafs")
         INPUT_DATA=$DATA_DIR/wafs.${wafs_grid}.full.grb2
         cp $INPUT_DATA $WORK_CTL/fort.9
         cp $INPUT_DATA $WORK_TEST/fort.9
         opt=1 ;;
      "ipxwafs2")
         INPUT_DATA=$DATA_DIR/wafs.${wafs_grid}.full.bitmap.grb2
         cp $INPUT_DATA $WORK_CTL/fort.9
         cp $INPUT_DATA $WORK_TEST/fort.9
         opt=2 ;;
      "ipxwafs3")
         INPUT_DATA=$DATA_DIR/wafs.${wafs_grid}.full.bitmap.grb2
         cp $INPUT_DATA $WORK_CTL/fort.9
         cp $INPUT_DATA $WORK_TEST/fort.9
         opt=3 ;;
      esac

      ctl_failed=0
      test_failed=0

      cd $WORK_TEST
      ./ipxwafs_test_${bytesize} ${opt} > log
      status=$?
      if ((status != 0))
      then
        echo "** TEST RUN FAILED **"
        mv log log.full2thin.${routine}.grid${wafs_grid}.${bytesize}byte.failed
        reg_test_failed=1
        test_failed=1
      fi

      cd $WORK_CTL
      ./ipxwafs_ctl_${bytesize} ${opt} > log
      status=$?
      if ((status != 0))
      then
        echo "** CONTROL RUN FAILED **"
        mv log log.full2thin.${routine}.grid${wafs_grid}.${bytesize}byte.failed
        reg_test_failed=1
        ctl_failed=1
      fi

      if ((ctl_failed == 0 && test_failed == 0));then
        cmp $WORK_CTL/output.grb2 $WORK_TEST/output.grb2
        status=$?
        if ((status != 0))
        then
          echo "** GRIB FILES NOT BIT IDENTICAL. REGRESSION TEST FAILED."
          mv $WORK_CTL/log $WORK_CTL/log.full2thin.${routine}.grid${wafs_grid}.${bytesize}byte.failed
          mv $WORK_CTL/output.grb2 $WORK_CTL/output.full2thin.${routine}.grid${wafs_grid}.${bytesize}byte.failed.grb2
          mv $WORK_TEST/log $WORK_TEST/log.full2thin.${routine}.grid${wafs_grid}.${bytesize}byte.failed
          mv $WORK_TEST/output.grb2 $WORK_TEST/output.full2thin.${routine}.grid${wafs_grid}.${bytesize}byte.failed.grb2
          reg_test_failed=1
        fi
      fi

      rm -f $WORK_CTL/output.grb2
      rm -f $WORK_CTL/log
      rm -f $WORK_CTL/fort.9
      rm -f $WORK_TEST/output.grb2
      rm -f $WORK_TEST/log
      rm -f $WORK_TEST/fort.9

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
