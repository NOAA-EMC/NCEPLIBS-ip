#!/bin/ksh

#------------------------------------------------------------------
# Sample script to run the ipolates regression test on
# zeus using 4 threads.  Modify path names as necessary.
#------------------------------------------------------------------

#@job_name=ipolates_4threads
#@output=ipolates_4threads.out
#@error=ipolates_4threads.out
#@class=dev
#@group=devonprod
#@account_no=GDAS-MTN
#@resources=ConsumableMemory(1000Mb)
#@job_type=parallel
#@task_affinity=cpu(4)
#@parallel_threads=4
#@wall_clock_limit=00:45:00
#@node_usage=shared
#@queue

# directory where regression tests reside
export REG_DIR=/global/save/wx20gg/gayno_ip_reg_tests/reg_tests

export WORK_DIR=/ptmp/$LOGNAME/regression
LOG_FILE=$WORK_DIR/ipolates_4threads.log

mkdir -p $WORK_DIR

$REG_DIR/ipolates/scripts/runall.ksh 4 >> $LOG_FILE

exit 0
