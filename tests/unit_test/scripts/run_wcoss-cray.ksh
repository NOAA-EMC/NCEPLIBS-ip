#!/bin/ksh
 
#-------------------------------------------------------
# Script to run the unit test on WCOSS-Cray
# compute nodes.
#
# Simply invoke this script on the command line
# with no arguments.
#
# Output is put in "unit_test.log"
#-------------------------------------------------------

set -x

bsub -oo unit_test.log -eo unit_test.log -q dev_shared -J ip_unit_test \
     -R rusage[mem=500] -P GFS-T2O -W 0:15 -cwd $(pwd) \
     "export OMP_NUM_THREADS=1; run_unit_test.ksh"

exit 0
