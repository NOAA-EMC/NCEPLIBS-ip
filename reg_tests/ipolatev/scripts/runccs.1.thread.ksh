#!/bin/ksh

#------------------------------------------------------------------
# Sample script to run the ipolatev regression test on
# ccs using 1 thread.  Modify path names as necessary.
#------------------------------------------------------------------

#@job_name=ipolatev_1thread
#@output=ipolatev_1thread.out
#@error=ipolatev_1thread.out
#@class=dev
#@group=devonprod
#@account_no=GDAS-MTN
#@resources=ConsumableMemory(1000Mb)
#@job_type=parallel
#@task_affinity=cpu(1)
#@parallel_threads=1
#@wall_clock_limit=00:45:00
#@node_usage=shared
#@queue

# directory where regression tests reside
export REG_DIR=/global/save/wx20gg/gayno_ip_reg_tests/reg_tests

export WORK_DIR=/ptmp/$LOGNAME/regression
LOG_FILE=$WORK_DIR/ipolatev_1thread.log

mkdir -p $WORK_DIR

$REG_DIR/ipolatev/scripts/runall.ksh 1 >> $LOG_FILE

exit 0
