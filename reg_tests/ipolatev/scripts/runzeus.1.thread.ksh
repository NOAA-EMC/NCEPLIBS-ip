#!/bin/ksh --login
#PBS -l nodes=1
#PBS -l walltime=1:00:00
# the account number. rm is regional model
#PBS -A rm
#PBS -N test
#PBS -o /scratch1/portfolios/NCEPDEV/da/save/George.Gayno/gayno_ip_reg_tests/reg_tests/ipolatev/scripts/log.1.thread
#PBS -e /scratch1/portfolios/NCEPDEV/da/save/George.Gayno/gayno_ip_reg_tests/reg_tests/ipolatev/scripts/log.1.thread

ulimit -s 1024000

export OMP_NUM_THREADS=1

cd /scratch1/portfolios/NCEPDEV/da/save/George.Gayno/gayno_ip_reg_tests/reg_tests/ipolatev/scripts

module load intel

./runall.ksh 1

exit 0
