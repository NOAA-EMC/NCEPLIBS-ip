#!/bin/ksh

#BSUB -eo ipolatev.log
#BSUB -oo ipolatev.log
#BSUB -q dev_shared
#BSUB -P GFS-T2O
#BSUB -a openmp -n 1
#BSUB -J ipolatev1
#BSUB -R affinity[core]
#BSUB -R rusage[mem=500]
#BSUB -R span[ptile=1]
#BSUB -W 0:20

set -x

export REG_DIR=${PWD}/../..

./runall.ksh 1

exit 0
