#!/bin/ksh

#------------------------------------------------------------------
# Sample script to run the ipolates regression test on
# ccs using 4 threads.  Modify path names as necessary.
#
# To run, type "llsubmit runccs.4.threads.ksh"
#------------------------------------------------------------------

#@job_name=ipolates_4threads
#@output=./log.4.threads
#@error=./log.4.threads
#@class=dev
#@group=dev
#@account_no=GDAS-MTN
#@resources=ConsumableMemory(1000Mb)
#@job_type=parallel
#@task_affinity=cpu(4)
#@parallel_threads=4
#@wall_clock_limit=00:45:00
#@node_usage=shared
#@queue

# directory where regression tests reside
export REG_DIR=../..

export WORK_DIR=/ptmp/$LOGNAME/regression
mkdir -p $WORK_DIR

./runall.ksh 4

exit 0
