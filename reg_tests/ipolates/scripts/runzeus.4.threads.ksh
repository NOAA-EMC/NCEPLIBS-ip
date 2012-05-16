#!/bin/ksh --login
#PBS -l nodes=1
#PBS -l walltime=0:30:00
# the account number. rm is regional model
#PBS -A rm
#PBS -N test
#PBS -o /scratch1/portfolios/NCEPDEV/da/save/George.Gayno/gayno_ip_reg_tests/reg_tests/ipolates/scripts/log.4.threads
#PBS -e /scratch1/portfolios/NCEPDEV/da/save/George.Gayno/gayno_ip_reg_tests/reg_tests/ipolates/scripts/log.4.threads

ulimit -s 1024000

export OMP_NUM_THREADS=4

cd /scratch1/portfolios/NCEPDEV/da/save/George.Gayno/gayno_ip_reg_tests/reg_tests/ipolates/scripts

module load intel

./runall.ksh 4

exit 0
