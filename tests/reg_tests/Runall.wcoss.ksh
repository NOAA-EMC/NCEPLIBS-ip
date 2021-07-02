#!/bin/ksh 

#------------------------------------------------------------------------
# Run the entire suite of ipolates2 (or ip2lib) regression tests on 
# the NCEP WCOSS phase 1/2 machine.
#
# See the README file for information on setting up and compiling
# the test suite.
#
# To run, type: "Runall.wcoss.ksh".  A series of "daisy-chained"
# job steps will be submitted.  To check the queue, type "bjobs".
#
# The run output is stored in $WORK_DIR.  Log output from the test suite 
# will be in "regression.log"  To monitor as the suite is running,
# do: grep ">>>" regression.log.  Once the suite is complete, a summary
# is placed in "summary.log" 
#------------------------------------------------------------------------

set -x

. /usrx/local/Modules/default/init/ksh
module load ics/15.0.6
module load lsf/9.1 

export REG_DIR=$(pwd)

export WORK_DIR="/stmpp1/${LOGNAME}/regression"
rm -fr $WORK_DIR
mkdir -p $WORK_DIR

LOG_FILE=${WORK_DIR}/regression.log
SUM_FILE=${WORK_DIR}/summary.log

bsub -e $LOG_FILE -o $LOG_FILE -q "dev_shared" -P "GFS-T2O" -J "gdswzd" -a openmp -n 1 \
     -R affinity[core] -R "rusage[mem=300]" -R span[ptile=1] -W 0:05 $REG_DIR/gdswzd/scripts/runall.ksh 

bsub -e $LOG_FILE -o $LOG_FILE -q "dev_shared" -P "GFS-T2O" \
     -J "ipxetas" -R affinity[core] -R "rusage[mem=100]" -W 0:05 -w 'ended(gdswzd)' $REG_DIR/ipxetas/scripts/runall.ksh 

bsub -e $LOG_FILE -o $LOG_FILE -q "dev_shared" -P "GFS-T2O" \
     -J "ipxwafs" -R affinity[core] -R "rusage[mem=100]" -W 0:05 -w 'ended(ipxetas)' $REG_DIR/ipxwafs/scripts/runall.ksh 

bsub -e $LOG_FILE -o $LOG_FILE -q "dev_shared" -P "GFS-T2O" -a openmp -n 1 \
     -J "ipolates1" -R affinity[core] -R "rusage[mem=500]" -R span[ptile=1] \
     -W 0:30 -w 'ended(ipxwafs)' $REG_DIR/ipolates/scripts/runall.ksh 1

bsub -e $LOG_FILE -o $LOG_FILE -q "dev_shared" -P "GFS-T2O" -a openmp -n 4 \
     -J "ipolates4" -R affinity[core] -R "rusage[mem=300]" -R span[ptile=4] \
     -W 0:30 -w 'ended(ipolates1)' $REG_DIR/ipolates/scripts/runall.ksh 4

bsub -e $LOG_FILE -o $LOG_FILE -q "dev_shared" -P "GFS-T2O" \
     -J "compares" -R affinity[core] -R "rusage[mem=100]" -W 0:10 -w 'ended(ipolates4)' $REG_DIR/ipolates/scripts/compare.ksh

bsub -e $LOG_FILE -o $LOG_FILE -q "dev_shared" -P "GFS-T2O" -a openmp -n 1 \
     -J "ipolatev1" -R affinity[core] -R "rusage[mem=500]" -R span[ptile=1] \
     -W 1:00 -w 'ended(compares)' $REG_DIR/ipolatev/scripts/runall.ksh 1

bsub -e $LOG_FILE -o $LOG_FILE -q "dev_shared" -P "GFS-T2O" -a openmp -n 4 \
     -J "ipolatev4" -R affinity[core] -R "rusage[mem=300]" -R span[ptile=4] \
     -W 1:00 -w 'ended(ipolatev1)' $REG_DIR/ipolatev/scripts/runall.ksh 4

bsub -e $LOG_FILE -o $LOG_FILE -q "dev_shared" -P "GFS-T2O" \
     -J "comparev" -R affinity[core] -R "rusage[mem=100]" -W 0:10 -w 'ended(ipolatev4)' $REG_DIR/ipolatev/scripts/compare.ksh

bsub -o $LOG_FILE -q "dev_shared" -P "GFS-T2O" -J "summary" \
     -R affinity[core] -R "rusage[mem=100]" -W 0:01 -w 'ended(comparev)' "grep '<<<' $LOG_FILE >> $SUM_FILE"

exit 0
