#!/bin/ksh --login
#
#-----------------------------------------------------------------------------
# Script to run the copygb regression test on Zeus.
#
# Invoke script by typing its name on the command line: "qsub run.zeus.ksh"
#-----------------------------------------------------------------------------
#
#PBS -l procs=1
#PBS -l mem=2500Mb
#PBS -l walltime=4:00:00
#PBS -A rm
#PBS -N iplib
#PBS -o ./regression.log
#PBS -e ./regression.log

set -x

export REG_DIR=$PBS_O_WORKDIR/../../

export WORK_DIR=/stmp/$LOGNAME/regression
mkdir -p $WORK_DIR

module load intel

ulimit -s 2048000

$REG_DIR/copygb/scripts/runall.ksh

exit 0
