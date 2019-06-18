#!/bin/ksh --login

#------------------------------------------------------------------
# Sample script to run the ipolates regression test on
# zeus using 1 thread.  Modify path names as necessary.
#
# To run, type "qsub runzeus.1.thread.ksh"
#------------------------------------------------------------------

#PBS -l nodes=1
#PBS -l walltime=0:30:00
# the account number. rm is regional model
#PBS -A rm
#PBS -N test
#PBS -o ./log.1.thread
#PBS -e ./log.1.thread

ulimit -s 2048000
ulimit -c unlimited

export OMP_NUM_THREADS=1

module load intel

export REG_DIR=/u/George.Gayno/save/ip_wcoss_port/reg_tests

export WORK_DIR=/u/George.Gayno/stmp/regression2
mkdir -p $WORK_DIR

${REG_DIR}/ipolates/scripts/runall.ksh 1

exit 0
