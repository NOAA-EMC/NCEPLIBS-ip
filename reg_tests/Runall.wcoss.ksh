#!/bin/ksh 

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
module load lsf/9.1 

export REG_DIR=$(pwd)

export WORK_DIR="/stmp/${LOGNAME}/regression"
rm -fr $WORK_DIR
mkdir -p $WORK_DIR

LOG_FILE=${WORK_DIR}/regression.log
SUM_FILE=${WORK_DIR}/summary.log

bsub -e $LOG_FILE -o $LOG_FILE -q "dev_shared" -P "GFS-MTN" \
     -J "gausslat" -R affinity[core] -R "rusage[mem=100]" -W 0:01 $REG_DIR/gausslat/scripts/runall.ksh 

bsub -e $LOG_FILE -o $LOG_FILE -q "dev_shared" -P "GFS-MTN" \
     -J "gcdist" -R affinity[core] -R "rusage[mem=100]" -W 0:01 -w 'ended(gausslat)' $REG_DIR/gcdist/scripts/runall.ksh 

bsub -e $LOG_FILE -o $LOG_FILE -q "dev_shared" -P "GFS-MTN" \
     -J "gdswiz" -R affinity[core] -R "rusage[mem=300]" -W 0:05 -w 'ended(gcdist)' $REG_DIR/gdswiz_wzd/scripts/runall.ksh 

bsub -e $LOG_FILE -o $LOG_FILE -q "dev_shared" -P "GFS-MTN" \
     -J "ipmerge2" -R affinity[core] -R "rusage[mem=100]" -W 0:02 -w 'ended(gdswiz)' $REG_DIR/ipmerge2/scripts/runall.ksh 

bsub -e $LOG_FILE -o $LOG_FILE -q "dev_shared" -P "GFS-MTN" \
     -J "ipsector" -R affinity[core] -R "rusage[mem=100]" -W 0:02 -w 'ended(ipmerge2)' $REG_DIR/ipsector/scripts/runall.ksh 

bsub -e $LOG_FILE -o $LOG_FILE -q "dev_shared" -P "GFS-MTN" \
     -J "ipxetas" -R affinity[core] -R "rusage[mem=100]" -W 0:02 -w 'ended(ipsector)' $REG_DIR/ipxetas/scripts/runall.ksh 

bsub -e $LOG_FILE -o $LOG_FILE -q "dev_shared" -P "GFS-MTN" \
     -J "ipxwafs" -R affinity[core] -R "rusage[mem=100]" -W 0:05 -w 'ended(ipxetas)' $REG_DIR/ipxwafs/scripts/runall.ksh 

bsub -e $LOG_FILE -o $LOG_FILE -q "dev_shared" -P "GFS-MTN" \
     -J "ipxwafs23" -R affinity[core] -R "rusage[mem=100]" -W 0:05 -w 'ended(ipxwafs)' $REG_DIR/ipxwafs2_3/scripts/runall.ksh 

bsub -e $LOG_FILE -o $LOG_FILE -q "dev_shared" -P "GFS-MTN" \
     -J "makgds" -R affinity[core] -R "rusage[mem=100]" -W 0:02 -w 'ended(ipxwafs23)' $REG_DIR/makgds/scripts/runall.ksh 

bsub -e $LOG_FILE -o $LOG_FILE -q "dev_shared" -P "GFS-MTN" -a openmp -n 1 \
     -J "ipolates1" -R affinity[core] -R "rusage[mem=500]" -R span[ptile=1] \
     -W 0:30 -w 'ended(makgds)' $REG_DIR/ipolates/scripts/runall.ksh 1

bsub -e $LOG_FILE -o $LOG_FILE -q "dev_shared" -P "GFS-MTN" -a openmp -n 4 \
     -J "ipolates4" -R affinity[core] -R "rusage[mem=300]" -R span[ptile=4] \
     -W 0:30 -w 'ended(ipolates1)' $REG_DIR/ipolates/scripts/runall.ksh 4

bsub -e $LOG_FILE -o $LOG_FILE -q "dev_shared" -P "GFS-MTN" \
     -J "compares" -R affinity[core] -R "rusage[mem=100]" -W 0:10 -w 'ended(ipolates4)' $REG_DIR/ipolates/scripts/compare.ksh

bsub -e $LOG_FILE -o $LOG_FILE -q "dev_shared" -P "GFS-MTN" -a openmp -n 1 \
     -J "ipolatev1" -R affinity[core] -R "rusage[mem=500]" -R span[ptile=1] \
     -W 1:00 -w 'ended(compares)' $REG_DIR/ipolatev/scripts/runall.ksh 1

bsub -e $LOG_FILE -o $LOG_FILE -q "dev_shared" -P "GFS-MTN" -a openmp -n 4 \
     -J "ipolatev4" -R affinity[core] -R "rusage[mem=300]" -R span[ptile=4] \
     -W 1:00 -w 'ended(ipolatev1)' $REG_DIR/ipolatev/scripts/runall.ksh 4

bsub -e $LOG_FILE -o $LOG_FILE -q "dev_shared" -P "GFS-MTN" \
     -J "comparev" -R affinity[core] -R "rusage[mem=100]" -W 0:10 -w 'ended(ipolatev4)' $REG_DIR/ipolatev/scripts/compare.ksh

bsub -o $LOG_FILE -q "dev_shared" -P "GFS-MTN" -J "summary" \
     -R affinity[core] -R "rusage[mem=100]" -W 0:01 -w 'ended(comparev)' "grep '<<<' $LOG_FILE >> $SUM_FILE"

exit 0
