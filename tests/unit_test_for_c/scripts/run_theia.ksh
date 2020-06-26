#!/bin/ksh --login

#------------------------------------------------------------
# Script to run the 'c' unit test on Theia compute nodes.
#
# To run, type: 'sbatch $script'
#
# Output is put in "unit_test.log"
#------------------------------------------------------------

<<<<<<< HEAD
#PBS -l procs=1
#PBS -l vmem=100M
#PBS -l walltime=0:01:00
#PBS -A glbss
#PBS -N ip_unit_test
#PBS -o unit_test.log
#PBS -e unit_test.log
=======
#SBATCH -p shared
#SBATCH --mem=100M
#SBATCH -t 0:01:00
#SBATCH -A fv3-cpu
#SBATCH -J ip_c_unit_test
#SBATCH -q debug
#SBATCH -o unit_test.log
#SBATCH -e unit_test.log
>>>>>>> 9930866fd8678c5bedd406fad303d3342d0c3f18

set -x

module purge
module load intel/15.6.233

rundir=${SLURM_SUBMIT_DIR}
cd $rundir

./run_unit_test.ksh

exit 0
