#!/bin/ksh --login

#------------------------------------------------------------------
# Sample script to run the ipolates regression test on
# zeus using 1 thread.  Modify path names as necessary.
#------------------------------------------------------------------

#PBS -l nodes=1
#PBS -l walltime=0:30:00
# the account number. rm is regional model
#PBS -A rm
#PBS -N test
#PBS -o /scratch1/portfolios/NCEPDEV/da/save/George.Gayno/gayno_ip_reg_tests/reg_tests/ipolates/scripts/log.1.thread
#PBS -e /scratch1/portfolios/NCEPDEV/da/save/George.Gayno/gayno_ip_reg_tests/reg_tests/ipolates/scripts/log.1.thread

ulimit -s 1024000

export OMP_NUM_THREADS=1

cd /scratch1/portfolios/NCEPDEV/da/save/George.Gayno/gayno_ip_reg_tests/reg_tests/ipolates/scripts

module load intel

./runall.ksh 1

exit 0
