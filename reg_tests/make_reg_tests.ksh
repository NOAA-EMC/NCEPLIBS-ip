#!/bin/ksh

#-----------------------------------------------------------------------------
# This script compiles all regression tests on WCOSS or Zeus using the 
# intel IFORT compiler.
#
# Before compiling, you must link the 'control' and 'test' IP libraries
# to the ./lib subdirectory, or copy them manually.  Normally, the
# 'control' library is the OPS version or the head of the trunk.  And
# the 'test' library is your branch version.  The naming convention is:
#
# libip_ctl_4.a  (4 byte integer/4 byte float, control)
# libip_ctl_8.a  (8 byte integer/8 byte float, control)
# libip_ctl_d.a  (4 byte integer/8 byte float, control)
# libip_test_4.a (4 byte integer/4 byte float, test)
# libip_test_8.a (8 byte integer/8 byte float, test)
# libip_test_d.a (4 byte integer/8 byte float, test)
#
# Invoke this script as follows: "make_reg_tests.ksh"
#-----------------------------------------------------------------------------

set -x

. ./config-setup/ifort.setup

MAKE="gmake"

# location and names of libraries are machine dependent

for WHICHIP in ctl test; do
  for PRECISION in 4 8 d; do

    case $PRECISION in
      d) PRECISION2=4 ;;
      *) PRECISION2=$PRECISION ;;
    esac

    ./configure --prefix=${PWD} --enable-promote=${PRECISION} \
      LIBS="-lip_${WHICHIP}_${PRECISION} -lsp_v${SP_LIB_V}_${PRECISION} -lbacio_v${BACIO_LIB_V}_${PRECISION2} -lw3nco_v${W3NCO_LIB_V}_${PRECISION}"
    if [ $? -ne 0 ]; then
      set +x
      echo "$0: Error configuring for ${PRECISION}-byte ${WHICHIP} version build." >&2
      exit 2
    fi

    $MAKE clean
    $MAKE
    if [ $? -ne 0 ]; then
      set +x
      echo "$0: Error building for ${PRECISION}-byte ${WHICHIP} version build." >&2
      exit 3
    fi

    $MAKE install suffix="_${WHICHIP}_${PRECISION}"
    if [ $? -ne 0 ]; then
      set +x
      echo "$0: Error installing for ${PRECISION}-byte ${WHICHIP} version build." >&2
     exit 4
    fi

    mv config.log config_${WHICHIP}_${PRECISION}.log

  done
done
