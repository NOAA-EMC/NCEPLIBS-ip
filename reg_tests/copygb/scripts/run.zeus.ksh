#!/bin/ksh --login
#
#-----------------------------------------------------------------------------
# Script to run the copygb regression test on Zeus.
#
# Invoke script by typing its name on the command line: "qsub run.zeus.ksh"
#-----------------------------------------------------------------------------
#
#PBS -l procs=1
#PBS -l vmem=2500M
#PBS -l walltime=4:00:00
#PBS -A gm
#PBS -N iplib
#PBS -o ./regression.log
#PBS -e ./regression.log

set -x

export REG_DIR=$PBS_O_WORKDIR/../../

export WORK_DIR=/scratch2/portfolios/NCEPDEV/stmp/$LOGNAME/regression
mkdir -p $WORK_DIR

module load intel

ulimit -s 2048000

export NWPROD=/scratch2/portfolios/NCEPDEV/global/save/glopara/svn/gfs/tags/REL-9.1.3/para/fix/fix_am

$REG_DIR/copygb/scripts/copygb.ksh

exit 0
