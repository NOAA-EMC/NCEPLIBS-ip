#!/bin/ksh --login

#----------------------------------------------------------------------------
# Run the entire suite of IPOLATES (or IPLIB) regression tests on Zeus.
#
# See the README file for information on setting up and compiling
# the test suite.
#
# Before, running set the $PROJECT_CODE to the project that will
# be charged when running the test suite.  To find out which
# projects you are authorized to use, type "account_params".
#
# To run, type:  "Runall.zeus.ksh". A series of "daisy-chained"
# job steps will be submitted.  To check the queue, type:
# "showq -n -v -u USERNAME"
#
# The run output is stored in $WORK_DIR.  Log output from the test suite
# will be in "regression.log"  To monitor as the suite is running,
# do: grep ">>>" regression.log.  Once the suite is complete, a summary
# is placed in "summary.log"
#----------------------------------------------------------------------------

# The project that will be charged when running these jobs.
PROJECT_CODE="rm"

# Location of the regression test directory.
export REG_DIR=$(pwd)

# Working directory.
export WORK_DIR=/scratch2/portfolios/NCEPDEV/stmp/$LOGNAME/regression
rm -fr $WORK_DIR
mkdir -p $WORK_DIR

# Output log files.
LOG_FILE=${WORK_DIR}/regression.log
SUM_FILE=${WORK_DIR}/summary.log

. /contrib/module/3.2.9/Modules/3.2.9/init/ksh
module load intel

GAUSSLAT=$(qsub -l nodes=1 -l walltime=0:01:00 -A $PROJECT_CODE -N iplib -o $LOG_FILE -e $LOG_FILE \
      -v REG_DIR,WORK_DIR $REG_DIR/gausslat/scripts/runall.ksh)

IPXETAS=$(qsub -l nodes=1 -l walltime=0:01:00 -A $PROJECT_CODE -N iplib -o $LOG_FILE -e $LOG_FILE \
      -v REG_DIR,WORK_DIR -W depend=afterok:$GAUSSLAT $REG_DIR/ipxetas/scripts/runall.ksh)

IPXWAFS=$(qsub -l nodes=1 -l walltime=0:02:00 -A $PROJECT_CODE -N iplib -o $LOG_FILE -e $LOG_FILE \
      -v REG_DIR,WORK_DIR -W depend=afterok:$IPXETAS $REG_DIR/ipxwafs/scripts/runall.ksh)

IPXWAFS2_3=$(qsub -l nodes=1 -l walltime=0:02:00 -A $PROJECT_CODE -N iplib -o $LOG_FILE -e $LOG_FILE \
      -v REG_DIR,WORK_DIR -W depend=afterok:$IPXWAFS $REG_DIR/ipxwafs2_3/scripts/runall.ksh)

MAKGDS=$(qsub -l nodes=1 -l walltime=0:02:00 -A $PROJECT_CODE -N iplib -o $LOG_FILE -e $LOG_FILE \
      -v REG_DIR,WORK_DIR -W depend=afterok:$IPXWAFS2_3 $REG_DIR/makgds/scripts/runall.ksh)

GDSWIZ=$(qsub -l nodes=1 -l walltime=0:30:00 -A $PROJECT_CODE -N iplib -o $LOG_FILE -e $LOG_FILE \
      -v REG_DIR,WORK_DIR -W depend=afterok:$MAKGDS $REG_DIR/gdswiz_wzd/scripts/runall.ksh)

export OMP_NUM_THREADS=1
IPOLATES_1=$(qsub -l nodes=1 -l walltime=1:00:00 -A $PROJECT_CODE -N iplib2 -o $LOG_FILE -e $LOG_FILE \
      -F "1" -v REG_DIR,WORK_DIR,OMP_NUM_THREADS -W depend=afterok:$GDSWIZ $REG_DIR/ipolates/scripts/runall.ksh)

export OMP_NUM_THREADS=4
IPOLATES_4=$(qsub -l nodes=1 -l walltime=1:00:00 -A $PROJECT_CODE -N iplib2 -o $LOG_FILE -e $LOG_FILE \
      -F "4" -W depend=afterok:$IPOLATES_1 \
      -v REG_DIR,WORK_DIR,OMP_NUM_THREADS $REG_DIR/ipolates/scripts/runall.ksh)

IPOLATES_CMP=$(qsub -l nodes=1 -l walltime=0:05:00 -A $PROJECT_CODE -N iplib2 -o $LOG_FILE -e $LOG_FILE \
      -v WORK_DIR -W depend=afterok:$IPOLATES_4 $REG_DIR/ipolates/scripts/compare.ksh)

export OMP_NUM_THREADS=1
IPOLATEV_1=$(qsub -l nodes=1 -l walltime=1:30:00 -A $PROJECT_CODE -N iplib2 -o $LOG_FILE -e $LOG_FILE \
      -F "1" -W depend=afterok:$IPOLATES_CMP \
      -v REG_DIR,WORK_DIR,OMP_NUM_THREADS $REG_DIR/ipolatev/scripts/runall.ksh)

export OMP_NUM_THREADS=4
IPOLATEV_4=$(qsub -l nodes=1 -l walltime=1:30:00 -A $PROJECT_CODE -N iplib2 -o $LOG_FILE -e $LOG_FILE \
      -F "4" -W depend=afterok:$IPOLATEV_1 \
      -v REG_DIR,WORK_DIR,OMP_NUM_THREADS $REG_DIR/ipolatev/scripts/runall.ksh)

IPOLATEV_CMP=$(qsub -l nodes=1 -l walltime=0:05:00 -A $PROJECT_CODE -N iplib2 -o $LOG_FILE -e $LOG_FILE \
      -v WORK_DIR -W depend=afterok:$IPOLATEV_4 $REG_DIR/ipolatev/scripts/compare.ksh)

SUMMARY=$(echo "grep '<<<' $LOG_FILE > $SUM_FILE" | qsub -l nodes=1 -l walltime=0:01:00 -A $PROJECT_CODE -N iplib -o $LOG_FILE -e $LOG_FILE \
      -v REG_DIR,WORK_DIR -W depend=afterok:$IPOLATEV_CMP)

exit 0
