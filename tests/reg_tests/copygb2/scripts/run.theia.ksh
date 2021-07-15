#!/bin/sh --login
#
#-----------------------------------------------------------------------------
# Script to run the copygb2 regression test on Theia.
#
# Invoke script by typing its name on the command line: "sbatch $script".
#-----------------------------------------------------------------------------
#
#SBATCH --ntasks=1
#SBATCH --mem=2500M
#SBATCH -t 1:00:00
#SBATCH -A fv3-cpu
#SBATCH -J iptest_copygb2
#SBATCH -q batch
#SBATCH -o ./regression.log
#SBATCH -e ./regression.log

set -x

module purge
module load intel

export REG_DIR=${SLURM_SUBMIT_DIR}/../../

export WORK_DIR=/scratch3/NCEPDEV/stmp1/$LOGNAME/regression
mkdir -p $WORK_DIR

export OMP_NUM_THREADS=1

$REG_DIR/copygb2/scripts/copygb2.ksh

exit 0
