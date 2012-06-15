#!/bin/ksh

#------------------------------------------------------------------
# Driver script to run the copygb regression test on CCS.
#
# Script is invoked by typing, "llsubmit run.ccs.ksh"
#------------------------------------------------------------------

#@job_name=copygb
#@output=regression.log
#@error=regression.log
#@class=dev
#@group=dev
#@account_no=GDAS-MTN
#@resources=ConsumableCpus(1)ConsumableMemory(2500Mb)
#@job_type=serial
#@wall_clock_limit=01:30:00
#@node_usage=shared
#@queue

# directory where regression tests reside.
export REG_DIR=$(pwd)/../..

# working directory.
export WORK_DIR=/ptmp/$LOGNAME/regression
mkdir -p $WORK_DIR

$REG_DIR/copygb/scripts/runall.ksh 

exit 0
