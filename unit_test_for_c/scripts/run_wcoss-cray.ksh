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
     -R rusage[mem=100] -P GFS-T2O -W 0:01 -cwd $(pwd) \
     "run_unit_test.ksh"

exit 0
