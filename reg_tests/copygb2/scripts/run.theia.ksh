#!/bin/ksh --login
#
#-----------------------------------------------------------------------------
# Script to run the copygb2 regression test on Theia.
#
# Invoke script by typing its name on the command line: "qsub run.theia.ksh"
#-----------------------------------------------------------------------------
#
#PBS -l procs=1
#PBS -l vmem=2500M
#PBS -l walltime=2:30:00
#PBS -A glbss
#PBS -N iptest_copygb2
#PBS -o ./regression.log
#PBS -e ./regression.log

set -x

module purge
module load intel

export REG_DIR=$PBS_O_WORKDIR/../../

export WORK_DIR=/scratch3/NCEPDEV/stmp1/$LOGNAME/regression
mkdir -p $WORK_DIR

export OMP_NUM_THREADS=1

$REG_DIR/copygb2/scripts/copygb2.ksh

exit 0
