#!/bin/ksh --login

#------------------------------------------------------------------
# Sample script to run the ipolates regression test on
# zeus using 4 threads.  Modify path names as necessary.
#
# To run, type "qsub runzeus.1.thread.ksh"
#------------------------------------------------------------------

#PBS -l nodes=1
#PBS -l walltime=0:30:00
# the account number. rm is regional model
#PBS -A rm
#PBS -N test
#PBS -o ./log.4.threads
#PBS -e ./log.4.threads

ulimit -s 2048000

export OMP_NUM_THREADS=4

module load intel

export REG_DIR=$PBS_O_WORKDIR/../../

export WORK_DIR=/stmp/$LOGNAME/regression
mkdir -p $WORK_DIR

$PBS_O_WORKDIR/runall.ksh 4

exit 0
