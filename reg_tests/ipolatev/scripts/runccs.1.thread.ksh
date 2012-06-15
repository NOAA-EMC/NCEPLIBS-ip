#!/bin/ksh

#------------------------------------------------------------------
# Sample script to run the ipolatev regression test on
# ccs using 1 thread.  Modify path names as necessary.
#
# To run, type "llsubmit runccs.1.thread.ksh"
#------------------------------------------------------------------

#@job_name=ipolatev_1thread
#@output=./log.1.thread
#@error=./log.1.thread
#@class=dev
#@group=dev
#@account_no=GDAS-MTN
#@resources=ConsumableMemory(1000Mb)
#@job_type=parallel
#@task_affinity=cpu(1)
#@parallel_threads=1
#@wall_clock_limit=01:15:00
#@node_usage=shared
#@queue

# directory where regression tests reside
export REG_DIR=../..

export WORK_DIR=/ptmp/$LOGNAME/regression
mkdir -p $WORK_DIR

./runall.ksh 1

exit 0
