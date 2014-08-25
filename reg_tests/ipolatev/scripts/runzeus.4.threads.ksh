#!/bin/ksh --login

#------------------------------------------------------------------
# Sample script to run the ipolatev regression test on
# zeus using 4 threads.  Modify path names as necessary.
#
# To run, type "qsub runzeus.4.threads.ksh"
#------------------------------------------------------------------

#PBS -l nodes=1:ppn=1
#PBS -l vmem=2000M
#PBS -l walltime=1:00:00
#PBS -A rm
#PBS -N ipolatev4
#PBS -o ./log.4.threads
#PBS -e ./log.4.threads

. /contrib/module/3.2.9/Modules/3.2.9/init/ksh
module load intel

export OMP_NUM_THREADS=4

export REG_DIR=$PBS_O_WORKDIR/../../

export WORK_DIR=/scratch2/portfolios/NCEPDEV/stmp/$LOGNAME/regression
mkdir -p $WORK_DIR

$PBS_O_WORKDIR/runall.ksh 4

exit 0
