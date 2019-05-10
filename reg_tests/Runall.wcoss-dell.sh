#!/bin/bash

#------------------------------------------------------------------------
# Run the entire suite of ipolates (or iplib) regression tests on 
# the NCEP WCOSS Phase 3 - Dell machine.
#
# See the README file for information on setting up and compiling
# the test suite.
#
# To run, type: "Runall.wcoss-dell.ksh".  A series of "daisy-chained"
# job steps will be submitted.  To check the queue, type "bjobs".
#
# The run output is stored in $WORK_DIR.  Log output from the test suite 
# will be in "regression.log"  To monitor as the suite is running,
# do: grep ">>>" regression.log.  Once the suite is complete, a summary
# is placed in "summary.log" 
#------------------------------------------------------------------------

set -x

module purge
module load EnvVars/1.0.2
module load ips/18.0.1.163
module load lsf/10.1

export REG_DIR=$(pwd)

export WORK_DIR="/gpfs/dell1/stmp/${LOGNAME}/regression"
rm -fr $WORK_DIR
mkdir -p $WORK_DIR

LOG_FILE=${WORK_DIR}/regression.log
SUM_FILE=${WORK_DIR}/summary.log

bsub -e $LOG_FILE -o $LOG_FILE -q "dev_shared" -P "GFS-T2O" \
     -J "gausslat" -R "affinity[core(1)]" -R "rusage[mem=100]" -W 0:01 -cwd $(pwd) $REG_DIR/gausslat/scripts/runall.ksh 

bsub -e $LOG_FILE -o $LOG_FILE -q "dev_shared" -P "GFS-T2O" \
     -J "gdswzd" -R "affinity[core(1)]" -R "rusage[mem=300]" -W 0:05 -w 'ended(gausslat)' -cwd $(pwd) $REG_DIR/gdswzd/scripts/runall.ksh 

bsub -e $LOG_FILE -o $LOG_FILE -q "dev_shared" -P "GFS-T2O" \
     -J "ipxwafs" -R "affinity[core(1)]" -R "rusage[mem=100]" -W 0:05 -w 'ended(gdswzd)' -cwd $(pwd) $REG_DIR/ipxwafs/scripts/runall.ksh 

bsub -e $LOG_FILE -o $LOG_FILE -q "dev_shared" -P "GFS-T2O" \
     -J "ipxwafs23" -R "affinity[core(1)]" -R "rusage[mem=100]" -W 0:05 -w 'ended(ipxwafs)' -cwd $(pwd) $REG_DIR/ipxwafs2_3/scripts/runall.ksh 

bsub -e $LOG_FILE -o $LOG_FILE -q "dev_shared" -P "GFS-T2O" \
     -J "makgds" -R "affinity[core(1)]" -R "rusage[mem=100]" -W 0:02 -w 'ended(ipxwafs23)' -cwd $(pwd) $REG_DIR/makgds/scripts/runall.ksh 

bsub -e $LOG_FILE -o $LOG_FILE -q "dev" -P "GFS-T2O" \
     -J "ipolates1" -R "rusage[mem=500]" -n 1 -R span[ptile=1] \
     -W 0:30 -w 'ended(makgds)' -cwd $(pwd) "export OMP_NUM_THREADS=1; $REG_DIR/ipolates/scripts/runall.ksh 1"

bsub -e $LOG_FILE -o $LOG_FILE -q "dev" -P "GFS-T2O"  \
     -J "ipolates4" -R "rusage[mem=300]" -n 4 -R span[ptile=4] \
     -W 0:30 -w 'ended(ipolates1)' -cwd $(pwd) "export OMP_NUM_THREADS=4; $REG_DIR/ipolates/scripts/runall.ksh 4"

bsub -e $LOG_FILE -o $LOG_FILE -q "dev_shared" -P "GFS-T2O" \
     -J "compares" -R "affinity[core(1)]" -R "rusage[mem=100]" -W 0:10 -w 'ended(ipolates4)' -cwd $(pwd) $REG_DIR/ipolates/scripts/compare.ksh

bsub -e $LOG_FILE -o $LOG_FILE -q "dev" -P "GFS-T2O" \
     -J "ipolatev1" -R "rusage[mem=500]" -n 1 -R span[ptile=1] \
     -W 1:00 -w 'ended(compares)' -cwd $(pwd) "export OMP_NUM_THREADS=1; $REG_DIR/ipolatev/scripts/runall.ksh 1"

bsub -e $LOG_FILE -o $LOG_FILE -q "dev" -P "GFS-T2O" \
     -J "ipolatev4" -R "rusage[mem=300]" -n 4 -R span[ptile=4] \
     -W 1:00 -w 'ended(ipolatev1)' -cwd $(pwd) "export OMP_NUM_THREADS=4; $REG_DIR/ipolatev/scripts/runall.ksh 4"

bsub -e $LOG_FILE -o $LOG_FILE -q "dev_shared" -P "GFS-T2O" \
     -J "comparev" -R "affinity[core(1)]" -R "rusage[mem=100]" -W 0:10 -w 'ended(ipolatev4)' -cwd $(pwd) $REG_DIR/ipolatev/scripts/compare.ksh

bsub -o $LOG_FILE -q "dev_shared" -P "GFS-T2O" -J "summary" \
     -R "affinity[core(1)]" -R "rusage[mem=100]" -W 0:01 -w 'ended(comparev)' -cwd $(pwd) "grep '<<<' $LOG_FILE >> $SUM_FILE"

exit 0
