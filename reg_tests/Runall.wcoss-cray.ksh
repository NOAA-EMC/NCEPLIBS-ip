#!/bin/ksh 

#------------------------------------------------------------------------
# Run the entire suite of ipolates (or iplib) regression tests on 
# the NCEP WCOSS-Cray machine.
#
# See the README file for information on setting up and compiling
# the test suite.
#
# To run, type: "Runall.wcoss-cray.ksh".  A series of "daisy-chained"
# job steps will be submitted.  To check the queue, type "bjobs".
#
# The run output is stored in $WORK_DIR.  Log output from the test suite 
# will be in "regression.log"  To monitor as the suite is running,
# do: grep ">>>" regression.log.  Once the suite is complete, a summary
# is placed in "summary.log" 
#------------------------------------------------------------------------

set -x

export REG_DIR=$(pwd)

export WORK_DIR="/gpfs/hps/stmp/${LOGNAME}/regression"
rm -fr $WORK_DIR
mkdir -p $WORK_DIR

LOG_FILE=${WORK_DIR}/regression.log
SUM_FILE=${WORK_DIR}/summary.log

bsub -e $LOG_FILE -o $LOG_FILE -q "dev_shared" -P "GFS-T2O" \
     -J "gausslat" -R "rusage[mem=100]" -W 0:01 -cwd $(pwd) $REG_DIR/gausslat/scripts/runall.ksh 

bsub -e $LOG_FILE -o $LOG_FILE -q "dev_shared" -P "GFS-T2O" \
     -J "gdswzd" -R "rusage[mem=300]" -W 0:05 -w 'ended(gausslat)' -cwd $(pwd) $REG_DIR/gdswzd/scripts/runall.ksh 

bsub -e $LOG_FILE -o $LOG_FILE -q "dev_shared" -P "GFS-T2O" \
     -J "ipxwafs" -R "rusage[mem=100]" -W 0:05 -w 'ended(gdswzd)' -cwd $(pwd) $REG_DIR/ipxwafs/scripts/runall.ksh 

bsub -e $LOG_FILE -o $LOG_FILE -q "dev_shared" -P "GFS-T2O" \
     -J "ipxwafs23" -R "rusage[mem=100]" -W 0:05 -w 'ended(ipxwafs)' -cwd $(pwd) $REG_DIR/ipxwafs2_3/scripts/runall.ksh 

bsub -e $LOG_FILE -o $LOG_FILE -q "dev_shared" -P "GFS-T2O" \
     -J "makgds" -R "rusage[mem=100]" -W 0:02 -w 'ended(ipxwafs23)' -cwd $(pwd) $REG_DIR/makgds/scripts/runall.ksh 

bsub -e $LOG_FILE -o $LOG_FILE -q "dev" -P "GFS-T2O" -J "ipolates1" \
     -M 500 -extsched 'CRAYLINUX[]' -W 0:30 -w 'ended(makgds)' -cwd $(pwd) \
      "export NODES=1; export OMP_NUM_THREADS=1; $REG_DIR/ipolates/scripts/runall.ksh 1"

bsub -e $LOG_FILE -o $LOG_FILE -q "dev" -P "GFS-T2O" -J "ipolates4" \
     -M 500 -extsched 'CRAYLINUX[]' -W 0:30 -w 'ended(ipolates1)' -cwd $(pwd) \
      "export NODES=1; export OMP_NUM_THREADS=4; $REG_DIR/ipolates/scripts/runall.ksh 4"

bsub -e $LOG_FILE -o $LOG_FILE -q "dev_shared" -P "GFS-T2O" \
     -J "compares" -R "rusage[mem=100]" -W 0:10 -w 'ended(ipolates4)' -cwd $(pwd) $REG_DIR/ipolates/scripts/compare.ksh

bsub -e $LOG_FILE -o $LOG_FILE -q "dev" -P "GFS-T2O" -J "ipolatev1" \
     -M 500 -extsched 'CRAYLINUX[]' -W 1:00 -w 'ended(compares)' -cwd $(pwd) \
      "export NODES=1; export OMP_NUM_THREADS=1; $REG_DIR/ipolatev/scripts/runall.ksh 1"

bsub -e $LOG_FILE -o $LOG_FILE -q "dev" -P "GFS-T2O" -J "ipolatev4" \
     -M 500 -extsched 'CRAYLINUX[]' -W 1:00 -w 'ended(ipolatev1)' -cwd $(pwd) \
      "export NODES=1; export OMP_NUM_THREADS=4; $REG_DIR/ipolatev/scripts/runall.ksh 4"

bsub -e $LOG_FILE -o $LOG_FILE -q "dev_shared" -P "GFS-T2O" \
     -J "comparev" -R "rusage[mem=100]" -W 0:10 -w 'ended(ipolatev4)' -cwd $(pwd) $REG_DIR/ipolatev/scripts/compare.ksh

bsub -o $LOG_FILE -q "dev_shared" -P "GFS-T2O" -J "summary" \
     -R "rusage[mem=100]" -W 0:01 -w 'ended(comparev)' -cwd $(pwd) "grep '<<<' $LOG_FILE >> $SUM_FILE"

exit 0
