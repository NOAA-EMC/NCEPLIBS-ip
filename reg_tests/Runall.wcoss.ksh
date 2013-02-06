#!/bin/ksh --login

#------------------------------------------------------------------------
# Run the entire series of iplib regression tests on the wcoss machine.
#
# To run, type: "Runall.wcoss.ksh"
#
# The log output, "regression.log" is stored in $WORK_DIR.
# A summary of the results will be placed in "summary.log" in $WORK_DIR.
#------------------------------------------------------------------------

set -x

. /usrx/local/Modules/3.2.9/init/ksh
module load ics/12.1
module load lsf/8.3

export REG_DIR=$(pwd)

export WORK_DIR="/stmp/${LOGNAME}/regression"
rm -fr $WORK_DIR
mkdir -p $WORK_DIR

LOG_FILE=${WORK_DIR}/regression.log
SUM_FILE=${WORK_DIR}/summary.log

bsub -e $LOG_FILE -o $LOG_FILE -q "dev" \
     -J "gausslat" -W 0:01 $REG_DIR/gausslat/scripts/runall.ksh 

bsub -e $LOG_FILE -o $LOG_FILE -q "dev" \
     -J "gcdist" -W 0:01 -w 'ended(gausslat)' $REG_DIR/gcdist/scripts/runall.ksh 

bsub -e $LOG_FILE -o $LOG_FILE -q "dev" \
     -J "gdswiz" -W 0:05 -w 'ended(gcdist)' $REG_DIR/gdswiz_wzd/scripts/runall.ksh 

bsub -e $LOG_FILE -o $LOG_FILE -q "dev" \
     -J "ipmerge2" -W 0:02 -w 'ended(gdswiz)' $REG_DIR/ipmerge2/scripts/runall.ksh 

bsub -e $LOG_FILE -o $LOG_FILE -q "dev" \
     -J "ipsector" -W 0:02 -w 'ended(ipmerge2)' $REG_DIR/ipsector/scripts/runall.ksh 

bsub -e $LOG_FILE -o $LOG_FILE -q "dev" \
     -J "ipxetas" -W 0:02 -w 'ended(ipsector)' $REG_DIR/ipxetas/scripts/runall.ksh 

bsub -e $LOG_FILE -o $LOG_FILE -q "dev" \
     -J "ipxwafs" -W 0:05 -w 'ended(ipxetas)' $REG_DIR/ipxwafs/scripts/runall.ksh 

bsub -e $LOG_FILE -o $LOG_FILE -q "dev" \
     -J "ipxwafs23" -W 0:05 -w 'ended(ipxwafs)' $REG_DIR/ipxwafs2_3/scripts/runall.ksh 

bsub -e $LOG_FILE -o $LOG_FILE -q "dev" \
     -J "makgds" -W 0:02 -w 'ended(ipxwafs23)' $REG_DIR/makgds/scripts/runall.ksh 

bsub -e $LOG_FILE -o $LOG_FILE -q "dev" -a openmp -n 1 \
     -J "ipolates1" -W 0:30 -w 'ended(makgds)' $REG_DIR/ipolates/scripts/runall.ksh 1

bsub -e $LOG_FILE -o $LOG_FILE -q "dev" -a openmp -n 4 \
     -J "ipolates4" -W 0:30 -w 'ended(ipolates1)' $REG_DIR/ipolates/scripts/runall.ksh 4

bsub -e $LOG_FILE -o $LOG_FILE -q "dev"  \
     -J "compares" -W 0:10 -w 'ended(ipolates4)' $REG_DIR/ipolates/scripts/compare.ksh

bsub -e $LOG_FILE -o $LOG_FILE -q "dev" -a openmp -n 1 \
     -J "ipolatev1" -W 1:00 -w 'ended(compares)' $REG_DIR/ipolatev/scripts/runall.ksh 1

bsub -e $LOG_FILE -o $LOG_FILE -q "dev" -a openmp -n 4 \
     -J "ipolatev4" -W 1:00 -w 'ended(ipolatev1)' $REG_DIR/ipolatev/scripts/runall.ksh 4

bsub -e $LOG_FILE -o $LOG_FILE -q "dev"  \
     -J "comparev" -W 0:10 -w 'ended(ipolatev4)' $REG_DIR/ipolatev/scripts/compare.ksh

bsub -o $LOG_FILE -q "dev" -J "summary" -W 0:01 -w 'ended(comparev)' "grep '<<<' $LOG_FILE >> $SUM_FILE"

exit 0
