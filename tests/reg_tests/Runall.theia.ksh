#!/bin/ksh --login

#----------------------------------------------------------------------------
# Run the entire suite of IPOLATES2 (or IP2LIB) regression tests on Theia.
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
# "squeue -u USERNAME"
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

GDSWZD=$(sbatch --parsable --ntasks=1 --mem=2000M -t 0:10:00 -A $PROJECT_CODE -J ip2test_gdswzd -o $LOG_FILE -e $LOG_FILE \
      -q batch --export REG_DIR,WORK_DIR,OMP_NUM_THREADS $REG_DIR/gdswzd/scripts/runall.ksh)

IPXETAS=$(sbatch --parsable --ntasks=1 --mem=500M -t 0:02:00 -A $PROJECT_CODE -J ip2test_ipxetas -o $LOG_FILE -e $LOG_FILE \
      --open-mode=append -q batch --export REG_DIR,WORK_DIR,OMP_NUM_THREADS -d afterok:$GDSWZD $REG_DIR/ipxetas/scripts/runall.ksh)

IPXWAFS=$(sbatch --parsable --ntasks=1 --mem=500M -t 0:02:00 -A $PROJECT_CODE -J ip2test_ipxwafs -o $LOG_FILE -e $LOG_FILE \
      --open-mode=append -q batch --export REG_DIR,WORK_DIR,OMP_NUM_THREADS -d afterok:$IPXETAS $REG_DIR/ipxwafs/scripts/runall.ksh)

IPOLATES_1=$(sbatch --parsable --ntasks=1 --mem=2000M -t 0:20:00 -A $PROJECT_CODE -J ip2test_ipolates1 -o $LOG_FILE -e $LOG_FILE \
      --open-mode=append -q batch --export REG_DIR,WORK_DIR,OMP_NUM_THREADS -d afterok:$IPXWAFS $REG_DIR/ipolates/scripts/runall.ksh "1")

export OMP_NUM_THREADS=4

IPOLATES_4=$(sbatch --parsable -N 1 -t 0:15:00 -A $PROJECT_CODE -J ip2test_ipolates4 -o $LOG_FILE -e $LOG_FILE \
      --open-mode=append -q batch -d afterok:$IPOLATES_1 \
      --export REG_DIR,WORK_DIR,OMP_NUM_THREADS $REG_DIR/ipolates/scripts/runall.ksh "4")

export OMP_NUM_THREADS=1

IPOLATES_CMP=$(sbatch --parsable --ntasks=1 --mem=200M -t 0:05:00 -A $PROJECT_CODE -J ip2test_ipolates_cmp -o $LOG_FILE -e $LOG_FILE \
      --open-mode=append -q batch --export WORK_DIR,OMP_NUM_THREADS -d afterok:$IPOLATES_4 $REG_DIR/ipolates/scripts/compare.ksh)

IPOLATEV_1=$(sbatch --parsable --ntasks=1 --mem=2000M -t 0:20:00 -A $PROJECT_CODE -J ip2test_ipolatev1 -o $LOG_FILE -e $LOG_FILE \
      --open-mode=append -q batch -d afterok:$IPOLATES_CMP \
      --export REG_DIR,WORK_DIR,OMP_NUM_THREADS $REG_DIR/ipolatev/scripts/runall.ksh "1")

export OMP_NUM_THREADS=4

IPOLATEV_4=$(sbatch --parsable -N 1 -t 0:15:00 -A $PROJECT_CODE -J ip2test_ipolatev4 -o $LOG_FILE -e $LOG_FILE \
      --open-mode=append -q batch -d afterok:$IPOLATEV_1 \
      --export REG_DIR,WORK_DIR,OMP_NUM_THREADS $REG_DIR/ipolatev/scripts/runall.ksh "4")

export OMP_NUM_THREADS=1

IPOLATEV_CMP=$(sbatch --parsable --ntasks=1 --mem=200M -t 0:05:00 -A $PROJECT_CODE -J ip2test_ipolatev_cmp -o $LOG_FILE -e $LOG_FILE \
      --open-mode=append -q batch --export WORK_DIR,OMP_NUM_THREADS -d afterok:$IPOLATEV_4 $REG_DIR/ipolatev/scripts/compare.ksh)

sbatch --ntasks=1 --mem=100M -t 0:01:00 -A $PROJECT_CODE -J ip2test_summary -o $LOG_FILE -e $LOG_FILE \
      --open-mode=append -q batch -d afterok:$IPOLATEV_CMP << EOF
#!/bin/sh
grep '<<<' $LOG_FILE > $SUM_FILE
EOF

exit 0
