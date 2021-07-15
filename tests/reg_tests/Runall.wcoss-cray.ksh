#!/bin/ksh 

#------------------------------------------------------------------------
# Run the entire suite of ipolates2 (or ip2lib) regression tests on 
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

export WORK_DIR="/gpfs/hps3/stmp/${LOGNAME}/regression"
rm -fr $WORK_DIR
mkdir -p $WORK_DIR

LOG_FILE=${WORK_DIR}/regression.log
SUM_FILE=${WORK_DIR}/summary.log

bsub -e $LOG_FILE -o $LOG_FILE -q "dev" -P "GFS-T2O" \
     -J "gdswzd" -M 500 -extsched 'CRAYLINUX[]' -W 0:10 \
     "export NODES=1; export OMP_NUM_THREADS=1; export machine=cray; $REG_DIR/gdswzd/scripts/runall.ksh"

bsub -e $LOG_FILE -o $LOG_FILE -q "debug" -P "GFS-T2O" \
     -J "ipxetas" -R "rusage[mem=100]" -W 0:05 -w 'ended(gdswzd)' $REG_DIR/ipxetas/scripts/runall.ksh 

bsub -e $LOG_FILE -o $LOG_FILE -q "debug" -P "GFS-T2O" \
     -J "ipxwafs" -R "rusage[mem=100]" -W 0:05 -w 'ended(ipxetas)' $REG_DIR/ipxwafs/scripts/runall.ksh 

bsub -e $LOG_FILE -o $LOG_FILE -q "dev" -P "GFS-T2O" -J "ipolates1" \
     -M 500 -extsched 'CRAYLINUX[]' -W 0:30 -w 'ended(ipxwafs)' \
      "export NODES=1; export OMP_NUM_THREADS=1; export machine=cray; $REG_DIR/ipolates/scripts/runall.ksh 1"

bsub -e $LOG_FILE -o $LOG_FILE -q "dev" -P "GFS-T2O" -J "ipolates4" \
     -M 500 -extsched 'CRAYLINUX[]' -W 0:30 -w 'ended(ipolates1)' \
      "export NODES=1; export OMP_NUM_THREADS=4; export machine=cray; $REG_DIR/ipolates/scripts/runall.ksh 4"

bsub -e $LOG_FILE -o $LOG_FILE -q "debug" -P "GFS-T2O" \
     -J "compares" -R "rusage[mem=100]" -W 0:10 -w 'ended(ipolates4)' $REG_DIR/ipolates/scripts/compare.ksh

bsub -e $LOG_FILE -o $LOG_FILE -q "dev" -P "GFS-T2O" -J "ipolatev1" \
     -M 500 -extsched 'CRAYLINUX[]' -W 1:00 -w 'ended(compares)' \
      "export NODES=1; export OMP_NUM_THREADS=1; export machine=cray; $REG_DIR/ipolatev/scripts/runall.ksh 1"

bsub -e $LOG_FILE -o $LOG_FILE -q "dev" -P "GFS-T2O" -J "ipolatev4" \
     -M 500 -extsched 'CRAYLINUX[]' -W 1:00 -w 'ended(ipolatev1)' \
      "export NODES=1; export OMP_NUM_THREADS=4; export machine=cray; $REG_DIR/ipolatev/scripts/runall.ksh 4"

bsub -e $LOG_FILE -o $LOG_FILE -q "debug" -P "GFS-T2O" \
     -J "comparev" -R "rusage[mem=100]" -W 0:10 -w 'ended(ipolatev4)' $REG_DIR/ipolatev/scripts/compare.ksh

bsub -o $LOG_FILE -q "debug" -P "GFS-T2O" -J "summary" \
     -R "rusage[mem=100]" -W 0:01 -w 'ended(comparev)' "grep '<<<' $LOG_FILE >> $SUM_FILE"

exit 0
