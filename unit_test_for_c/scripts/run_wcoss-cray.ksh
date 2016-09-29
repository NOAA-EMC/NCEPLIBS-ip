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

#BSUB -oo unit_test.log
#BSUB -eo unit_test.log
#BSUB -q dev_shared
#BSUB -J ip_unit_test
#BSUB -R rusage[mem=100]
#BSUB -P GFS-T2O
#BSUB -W 0:01

set -x

bsub -oo unit_test.log -eo unit_test.log -q dev_shared -J ip_unit_test \
     -R rusage[mem=100] -P GFS-T2O -W 0:01 -cwd $(pwd) \
     "run_unit_test.ksh"

exit 0
