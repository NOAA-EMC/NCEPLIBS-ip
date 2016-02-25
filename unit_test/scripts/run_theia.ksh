#!/bin/ksh --login

#------------------------------------------------------------
# To run, type 'qsub $script'
#------------------------------------------------------------

#PBS -l procs=1
#PBS -l vmem=500M
#PBS -l walltime=0:15:00
#PBS -A glbss
#PBS -N ip_unit_test
#PBS -o unit_test.log
#PBS -e unit_test.log

set -x

export OMP_NUM_THREADS=1

rundir=$PBS_O_WORKDIR
cd $rundir

run_unit_test.ksh

exit 0
