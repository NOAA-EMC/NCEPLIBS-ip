#!/bin/ksh --login

#----------------------------------------------------------------------------
# Run the entire suite of IPOLATES (or IPLIB) regression tests on Theia.
#
# See the README file for information on setting up and compiling
# the test suite.
#
# Before, running set the $PROJECT_CODE to the project that will
# be charged when running the test suite.  To find out which
# projects you are authorized to use, type "account_params".
#
# To run, type:  "Runall.theia.ksh". A series of "daisy-chained"
# job steps will be submitted.  To check the queue, type:
# "showq -n -v -u USERNAME"
#
# The run output is stored in $WORK_DIR.  Log output from the test suite
# will be in "regression.log"  To monitor as the suite is running,
# do: grep ">>>" regression.log.  Once the suite is complete, a summary
# is placed in "summary.log"
#----------------------------------------------------------------------------

# The project that will be charged when running these jobs.
PROJECT_CODE=${PROJECT_CODE:-"fv3-cpu"}

# Location of the regression test directory.
export REG_DIR=$(pwd)

# Working directory.
export WORK_DIR=/scratch3/NCEPDEV/stmp1/$LOGNAME/regression
rm -fr $WORK_DIR
mkdir -p $WORK_DIR

# Output log files.
LOG_FILE=${WORK_DIR}/regression.log
SUM_FILE=${WORK_DIR}/summary.log

module purge
module load intel

export OMP_NUM_THREADS=1

GAUSSLAT=$(qsub -l procs=1 -l vmem=500M -l walltime=0:01:00 -A $PROJECT_CODE -N iptest_gausslat -o $LOG_FILE -e $LOG_FILE \
      -v REG_DIR,WORK_DIR,OMP_NUM_THREADS $REG_DIR/gausslat/scripts/runall.ksh)

IPXWAFS=$(qsub -l procs=1 -l vmem=500M -l walltime=0:02:00 -A $PROJECT_CODE -N iptest_ipxwafs -o $LOG_FILE -e $LOG_FILE \
      -v REG_DIR,WORK_DIR,OMP_NUM_THREADS -W depend=afterok:$GAUSSLAT $REG_DIR/ipxwafs/scripts/runall.ksh)

IPXWAFS2_3=$(qsub -l procs=1 -l vmem=500M -l walltime=0:02:00 -A $PROJECT_CODE -N iptest_ipxwafs2 -o $LOG_FILE -e $LOG_FILE \
      -v REG_DIR,WORK_DIR,OMP_NUM_THREADS -W depend=afterok:$IPXWAFS $REG_DIR/ipxwafs2_3/scripts/runall.ksh)

MAKGDS=$(qsub -l procs=1 -l vmem=500M -l walltime=0:02:00 -A $PROJECT_CODE -N iptest_makgds -o $LOG_FILE -e $LOG_FILE \
      -v REG_DIR,WORK_DIR,OMP_NUM_THREADS -W depend=afterok:$IPXWAFS2_3 $REG_DIR/makgds/scripts/runall.ksh)

GDSWZD=$(qsub -l procs=1 -l vmem=2000M -l walltime=0:10:00 -A $PROJECT_CODE -N iptest_gdswzd -o $LOG_FILE -e $LOG_FILE \
      -v REG_DIR,WORK_DIR,OMP_NUM_THREADS -W depend=afterok:$MAKGDS $REG_DIR/gdswzd/scripts/runall.ksh)

IPOLATES_1=$(qsub -l procs=1 -l vmem=2000M -l walltime=0:30:00 -A $PROJECT_CODE -N iptest_ipolates1 -o $LOG_FILE -e $LOG_FILE \
      -F "1" -v REG_DIR,WORK_DIR,OMP_NUM_THREADS -W depend=afterok:$GDSWZD $REG_DIR/ipolates/scripts/runall.ksh)

export OMP_NUM_THREADS=4

IPOLATES_4=$(qsub -l nodes=1:ppn=24 -l walltime=0:30:00 -A $PROJECT_CODE -N iptest_ipolates4 -o $LOG_FILE -e $LOG_FILE \
      -F "4" -W depend=afterok:$IPOLATES_1 \
      -v REG_DIR,WORK_DIR,OMP_NUM_THREADS $REG_DIR/ipolates/scripts/runall.ksh)

export OMP_NUM_THREADS=1

IPOLATES_CMP=$(qsub -l procs=1 -l vmem=2000M -l walltime=0:05:00 -A $PROJECT_CODE -N iptest_ipolates_cmp -o $LOG_FILE -e $LOG_FILE \
      -v WORK_DIR,OMP_NUM_THREADS -W depend=afterok:$IPOLATES_4 $REG_DIR/ipolates/scripts/compare.ksh)

IPOLATEV_1=$(qsub -l procs=1 -l vmem=2000M -l walltime=0:45:00 -A $PROJECT_CODE -N iptest_ipolatev1 -o $LOG_FILE -e $LOG_FILE \
      -F "1" -W depend=afterok:$IPOLATES_CMP \
      -v REG_DIR,WORK_DIR,OMP_NUM_THREADS $REG_DIR/ipolatev/scripts/runall.ksh)

export OMP_NUM_THREADS=4

IPOLATEV_4=$(qsub -l nodes=1:ppn=24 -l walltime=0:30:00 -A $PROJECT_CODE -N iptest_ipolatev4 -o $LOG_FILE -e $LOG_FILE \
      -F "4" -W depend=afterok:$IPOLATEV_1 \
      -v REG_DIR,WORK_DIR,OMP_NUM_THREADS $REG_DIR/ipolatev/scripts/runall.ksh)

export OMP_NUM_THREADS=1

IPOLATEV_CMP=$(qsub -l procs=1 -l vmem=2000M -l walltime=0:05:00 -A $PROJECT_CODE -N iptest_ipolatev_cmp -o $LOG_FILE -e $LOG_FILE \
      -v WORK_DIR,OMP_NUM_THREADS -W depend=afterok:$IPOLATEV_4 $REG_DIR/ipolatev/scripts/compare.ksh)

SUMMARY=$(echo "grep '<<<' $LOG_FILE > $SUM_FILE" | qsub -l procs=1 -l vmem=500M -l walltime=0:01:00 -A $PROJECT_CODE -N iptest_summary \
      -o $LOG_FILE -e $LOG_FILE  -v REG_DIR,WORK_DIR,OMP_NUM_THREADS -W depend=afterok:$IPOLATEV_CMP)

exit 0
