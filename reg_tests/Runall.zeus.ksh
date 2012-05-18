#!/bin/ksh --login

set -x

export REG_DIR=$(pwd)

typeset -L4 machine
machine=$(dnsdomainname)

if [[ $machine == "zeus" ]];then    # zeus
  export WORK_DIR=${WORK_DIR:-/scratch2/portfolios/NCEPDEV/stmp/$LOGNAME/regression}
else   # cirrus/stratus
  export WORK_DIR=${WORK_DIR:-/stmp/$LOGNAME/regression}
fi

rm -fr $WORK_DIR
mkdir -p $WORK_DIR

LOG_FILE=${WORK_DIR}/regression.log

ulimit -s 1024000

module load intel

GAUSSLAT=$(qsub -l nodes=1 -l walltime=0:01:00 -A rm -N iplib -o $LOG_FILE -e $LOG_FILE \
      -v REG_DIR,WORK_DIR $REG_DIR/gausslat/scripts/runall.ksh)

GCDIST=$(qsub -l nodes=1 -l walltime=0:01:00 -A rm -N iplib -o $LOG_FILE -e $LOG_FILE \
      -v REG_DIR,WORK_DIR -W depend=afterok:$GAUSSLAT $REG_DIR/gcdist/scripts/runall.ksh)

IMERGE2=$(qsub -l nodes=1 -l walltime=0:01:00 -A rm -N iplib -o $LOG_FILE -e $LOG_FILE \
      -v REG_DIR,WORK_DIR -W depend=afterok:$GCDIST $REG_DIR/ipmerge2/scripts/runall.ksh)

IPSECTOR=$(qsub -l nodes=1 -l walltime=0:01:00 -A rm -N iplib -o $LOG_FILE -e $LOG_FILE \
      -v REG_DIR,WORK_DIR -W depend=afterok:$IMERGE2 $REG_DIR/ipsector/scripts/runall.ksh)

IPXETAS=$(qsub -l nodes=1 -l walltime=0:01:00 -A rm -N iplib -o $LOG_FILE -e $LOG_FILE \
      -v REG_DIR,WORK_DIR -W depend=afterok:$IPSECTOR $REG_DIR/ipxetas/scripts/runall.ksh)

IPXWAFS=$(qsub -l nodes=1 -l walltime=0:02:00 -A rm -N iplib -o $LOG_FILE -e $LOG_FILE \
      -v REG_DIR,WORK_DIR -W depend=afterok:$IPXETAS $REG_DIR/ipxwafs/scripts/runall.ksh)

IPXWAFS2_3=$(qsub -l nodes=1 -l walltime=0:02:00 -A rm -N iplib -o $LOG_FILE -e $LOG_FILE \
      -v REG_DIR,WORK_DIR -W depend=afterok:$IPXWAFS $REG_DIR/ipxwafs2_3/scripts/runall.ksh)

MAKGDS=$(qsub -l nodes=1 -l walltime=0:02:00 -A rm -N iplib -o $LOG_FILE -e $LOG_FILE \
      -v REG_DIR,WORK_DIR -W depend=afterok:$IPXWAFS2_3 $REG_DIR/makgds/scripts/runall.ksh)

GDSWIZ=$(qsub -l nodes=1 -l walltime=0:02:00 -A rm -N iplib -o $LOG_FILE -e $LOG_FILE \
      -v REG_DIR,WORK_DIR -W depend=afterok:$MAKGDS $REG_DIR/gdswiz_wzd/scripts/runall.ksh)

export OMP_NUM_THREADS=1
IPOLATES_1=$(qsub -l nodes=1 -l walltime=0:02:00 -A rm -N iplib2 -o $LOG_FILE -e $LOG_FILE \
      -F "1" -v REG_DIR,WORK_DIR,OMP_NUM_THREADS -W depend=afterok:$GDSWIZ $REG_DIR/ipolates/scripts/runall.ksh)

export OMP_NUM_THREADS=4
IPOLATES_4=$(qsub -l nodes=1 -l walltime=0:02:00 -A rm -N iplib2 -o $LOG_FILE -e $LOG_FILE \
      -F "4" -W depend=afterok:$IPOLATES_1 \
      -v REG_DIR,WORK_DIR,OMP_NUM_THREADS $REG_DIR/ipolates/scripts/runall.ksh)

IPOLATES_CMP=$(qsub -l nodes=1 -l walltime=0:05:00 -A rm -N iplib2 -o $LOG_FILE -e $LOG_FILE \
      -v WORK_DIR -W depend=afterok:$IPOLATES_4 $REG_DIR/ipolates/scripts/compare.ksh)

export OMP_NUM_THREADS=1
IPOLATEV_1=$(qsub -l nodes=1 -l walltime=0:30:00 -A rm -N iplib2 -o $LOG_FILE -e $LOG_FILE \
      -F "1" -W depend=afterok:$IPOLATES_CMP \
      -v REG_DIR,WORK_DIR,OMP_NUM_THREADS $REG_DIR/ipolatev/scripts/runall.ksh)

export OMP_NUM_THREADS=4
IPOLATEV_4=$(qsub -l nodes=1 -l walltime=0:30:00 -A rm -N iplib2 -o $LOG_FILE -e $LOG_FILE \
      -F "4" -W depend=afterok:$IPOLATEV_1 \
      -v REG_DIR,WORK_DIR,OMP_NUM_THREADS $REG_DIR/ipolatev/scripts/runall.ksh)

IPOLATEV_CMP=$(qsub -l nodes=1 -l walltime=0:05:00 -A rm -N iplib2 -o $LOG_FILE -e $LOG_FILE \
      -v WORK_DIR -W depend=afterok:$IPOLATEV_4 $REG_DIR/ipolatev/scripts/compare.ksh)

exit 0
