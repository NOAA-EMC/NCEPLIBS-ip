#!/bin/ksh --login

#------------------------------------------------------------------
# Sample script to run the ipolates regression test on
# zeus using 1 thread.  Modify path names as necessary.
#
# To run, type "qsub runzeus.1.thread.ksh"
#------------------------------------------------------------------

#PBS -l procs=1
#PBS -l vmem=2000M
#PBS -l walltime=0:30:00
#PBS -A rm
#PBS -N ipolates_1
#PBS -o ./log.1.thread
#PBS -e ./log.1.thread

. /contrib/module/3.2.9/Modules/3.2.9/init/ksh
module load intel

export OMP_NUM_THREADS=1

export REG_DIR=$PBS_O_WORKDIR/../../

export WORK_DIR=/scratch2/portfolios/NCEPDEV/stmp/$LOGNAME/regression
mkdir -p $WORK_DIR

$PBS_O_WORKDIR/runall.ksh 1

exit 0
