#!/bin/ksh --login
#
#-----------------------------------------------------------------------------
# Script to run the copygb regression test on Theia.
#
# Invoke script by typing its name on the command line: "qsub run.theia.ksh"
#-----------------------------------------------------------------------------
#
#PBS -l procs=1
#PBS -l vmem=2500M
#PBS -l walltime=2:30:00
#PBS -A fv3-cpu
#PBS -N iptest_copygb
#PBS -o ./regression.log
#PBS -e ./regression.log

set -x

module purge
module load intel

export REG_DIR=$PBS_O_WORKDIR/../../

export WORK_DIR=/scratch3/NCEPDEV/stmp1/$LOGNAME/regression
mkdir -p $WORK_DIR

export NWPROD=/scratch4/NCEPDEV/global/save/glopara/svn/gfs/tags/gfs.v13.0.2/para/fix/fix_am

export OMP_NUM_THREADS=1

$REG_DIR/copygb/scripts/copygb.ksh

exit 0
