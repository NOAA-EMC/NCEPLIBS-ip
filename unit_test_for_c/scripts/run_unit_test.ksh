#!/bin/ksh

#---------------------------------------------------------------------------------
# Driver script to run the IPOLATES (iplib) 'C' wrapper unit test.
#
# $Id$
#
# This script calls a 'C' program which calls iplib routine 'gdswzd' to
# compute the grid specs for a rotated lat/lon 'B'-grid.  The output is 
# sent to standard output.
#
# To run this script interactively, type "run_unit_test.sh".
# On WCOSS Phase 1/2, this script may be submitted to the compute nodes
# using "./run_wcoss.lsf".  On Theia, use "./run_theia.ksh"
# On WCOSS-Cray, use "./run_wcoss-cray.lsf".
#
# The source code for the 'C' program is located in ../sorc.
# There are separate versions for the single ("4"), mixed ("d"), and
# double ("8") versions of iplib.  After compilation, the executables
# are in ../exec.
#---------------------------------------------------------------------------------

#set -x

for  precision in "4" "d" "8"  # test all three precision versions of ipolates
do

  echo
  echo "********************************************************"
  echo "*** TEST $precision BYTE VERSION OF IPOLATES LIBRARY ***"
  echo "********************************************************"
  echo

  EXEC="../exec/test_gdswzd_${precision}.exe"

  $EXEC

done

echo
echo "**************"
echo "**** DONE ****"
echo "**************"

exit 0
