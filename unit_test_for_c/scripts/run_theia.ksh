#!/bin/ksh --login

#------------------------------------------------------------
# Script to run the 'c' unit test on Theia compute nodes.
#
# To run, type: 'qsub run_theia.ksh'
#
# Output is put in "unit_test.log"
#------------------------------------------------------------

#PBS -l procs=1
#PBS -l vmem=100M
#PBS -l walltime=0:01:00
#PBS -A fv3-cpu
#PBS -N ip_unit_test
#PBS -o unit_test.log
#PBS -e unit_test.log

set -x

rundir=$PBS_O_WORKDIR
cd $rundir

run_unit_test.ksh

exit 0
