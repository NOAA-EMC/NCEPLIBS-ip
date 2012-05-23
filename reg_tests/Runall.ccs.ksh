#!/bin/ksh

#------------------------------------------------------------------------
# Run the entire series of iplib regression tests on the ibm ccs machine.
#
# To run, type: "llsubmit Runall.ccs.ksh"
#
# The log output, "regression.log" is stored in $WORK_DIR.  
#------------------------------------------------------------------------

set -x

#@job_name=regression_driver
#@output=regress.out
#@error=regress.out
#@job_type=serial
#@class=dev
#@group=devonprod
#@account_no=GDAS-MTN

#@step_name=gausslat
#@job_type=serial
#@wall_clock_limit=00:01:00
#@node_usage=shared
#@resources=ConsumableCpus(1)ConsumableMemory(100Mb)
#@queue

#@step_name=ipolates_1thread
#@resources=ConsumableMemory(1000Mb)
#@job_type=parallel
#@task_affinity=cpu(1)
#@parallel_threads=1
#@wall_clock_limit=00:45:00
#@node_usage=shared
#@dependency=(gausslat == 0)
#@queue

#@step_name=ipolates_4thread
#@resources=ConsumableMemory(1000Mb)
#@job_type=parallel
#@task_affinity=cpu(4)
#@parallel_threads=4
#@wall_clock_limit=00:30:00
#@node_usage=shared
#@dependency=(ipolates_1thread == 0)
#@queue

#@step_name=ipolates_compare
#@resources=ConsumableMemory(1000Mb)
#@job_type=serial
#@task_affinity=cpu(1)
#@parallel_threads=1
#@wall_clock_limit=00:03:00
#@node_usage=shared
#@dependency=(ipolates_4thread == 0)
#@queue

#@step_name=ipolatev_1thread
#@resources=ConsumableMemory(1000Mb)
#@job_type=parallel
#@task_affinity=cpu(1)
#@parallel_threads=1
#@wall_clock_limit=00:30:00
#@node_usage=shared
#@dependency=(ipolates_compare == 0)
#@queue

#@step_name=ipolatev_4thread
#@resources=ConsumableMemory(1000Mb)
#@job_type=parallel
#@task_affinity=cpu(4)
#@parallel_threads=4
#@wall_clock_limit=00:20:00
#@node_usage=shared
#@dependency=(ipolatev_1thread == 0)
#@queue

#@step_name=ipolatev_compare
#@resources=ConsumableMemory(1000Mb)
#@job_type=serial
#@task_affinity=cpu(1)
#@parallel_threads=1
#@wall_clock_limit=00:03:00
#@node_usage=shared
#@dependency=(ipolatev_4thread == 0)
#@queue

#@step_name=gcdist
#@resources=ConsumableMemory(1000Mb)
#@job_type=serial
#@task_affinity=cpu(1)
#@parallel_threads=1
#@wall_clock_limit=00:03:00
#@node_usage=shared
#@dependency=(ipolatev_compare == 0)
#@queue

#@step_name=ipmerge2
#@resources=ConsumableMemory(1000Mb)
#@job_type=serial
#@task_affinity=cpu(1)
#@parallel_threads=1
#@wall_clock_limit=00:03:00
#@node_usage=shared
#@dependency=(gcdist == 0)
#@queue

#@step_name=ipsector
#@resources=ConsumableMemory(1000Mb)
#@job_type=serial
#@task_affinity=cpu(1)
#@parallel_threads=1
#@wall_clock_limit=00:03:00
#@node_usage=shared
#@dependency=(ipmerge2 == 0)
#@queue

#@step_name=ipxetas
#@resources=ConsumableMemory(1000Mb)
#@job_type=serial
#@task_affinity=cpu(1)
#@parallel_threads=1
#@wall_clock_limit=00:03:00
#@node_usage=shared
#@dependency=(ipsector == 0)
#@queue

#@step_name=ipxwafs
#@resources=ConsumableMemory(1000Mb)
#@job_type=serial
#@task_affinity=cpu(1)
#@parallel_threads=1
#@wall_clock_limit=00:03:00
#@node_usage=shared
#@dependency=(ipxetas == 0)
#@queue

#@step_name=ipxwafs2_3
#@resources=ConsumableMemory(1000Mb)
#@job_type=serial
#@task_affinity=cpu(1)
#@parallel_threads=1
#@wall_clock_limit=00:03:00
#@node_usage=shared
#@dependency=(ipxwafs == 0)
#@queue

#@step_name=gdswiz
#@resources=ConsumableMemory(1000Mb)
#@job_type=serial
#@task_affinity=cpu(1)
#@parallel_threads=1
#@wall_clock_limit=00:15:00
#@node_usage=shared
#@dependency=(ipxwafs2_3 == 0)
#@queue

#@step_name=makgds
#@resources=ConsumableMemory(1000Mb)
#@job_type=serial
#@task_affinity=cpu(1)
#@parallel_threads=1
#@wall_clock_limit=00:03:00
#@node_usage=shared
#@dependency=(gdswiz == 0)
#@queue

# directory where regression tests reside
export REG_DIR=$(pwd)

export WORK_DIR=/ptmp/$LOGNAME/regression
LOG_FILE=$WORK_DIR/regression.log

case $LOADL_STEP_NAME in
gausslat)
  rm -fr $WORK_DIR
  mkdir -p $WORK_DIR
  $REG_DIR/gausslat/scripts/runall.ksh > $LOG_FILE;;
gcdist)
  $REG_DIR/gcdist/scripts/runall.ksh >> $LOG_FILE;;
ipmerge2)
  $REG_DIR/ipmerge2/scripts/runall.ksh >> $LOG_FILE;;
ipsector)
  $REG_DIR/ipsector/scripts/runall.ksh >> $LOG_FILE;;
ipxetas)
  $REG_DIR/ipxetas/scripts/runall.ksh >> $LOG_FILE;;
ipxwafs)
  $REG_DIR/ipxwafs/scripts/runall.ksh >> $LOG_FILE;;
ipxwafs2_2)
  $REG_DIR/ipxwafs2_2/scripts/runall.ksh >> $LOG_FILE;;
gdswiz)
  $REG_DIR/gdswiz_wzd/scripts/runall.ksh >> $LOG_FILE;;
ipolates_1thread)
  $REG_DIR/ipolates/scripts/runall.ksh 1 >> $LOG_FILE;;
ipolates_4thread)
  $REG_DIR/ipolates/scripts/runall.ksh 4 >> $LOG_FILE;;
ipolates_compare)
  $REG_DIR/ipolates/scripts/compare.ksh >> $LOG_FILE;;
ipolatev_1thread)
  $REG_DIR/ipolatev/scripts/runall.ksh 1 >> $LOG_FILE;;
ipolatev_4thread)
  $REG_DIR/ipolatev/scripts/runall.ksh 4 >> $LOG_FILE;;
ipolatev_compare)
  $REG_DIR/ipolatev/scripts/compare.ksh >> $LOG_FILE;;
makgds)
  $REG_DIR/makgds/scripts/runall.ksh >> $LOG_FILE
  echo "ALL REGRESSION TESTS COMPLETED" >> $LOG_FILE;;
esac

exit 0
