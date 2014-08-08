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

#. /contrib/module/3.2.9/Modules/3.2.9/init/ksh
#module use -a /contrib/nceplibs/Modules/modulefiles
#module load sp

. ./config-setup/ifort.setup

BACIO_LIB4=${BACIO_LIB4:?}
SP_LIB4=${SP_LIB4:?}
W3NCO_LIB4=${W3NCO_LIB4:?}

BACIO_LIB8=${BACIO_LIB8:?}
SP_LIB8=${SP_LIB8:?}
W3NCO_LIB8=${W3NCO_LIB8:?}

SP_LIBd=${SP_LIBd:?}
W3NCO_LIBd=${W3NCO_LIBd:?}

MAKE="gmake"

# location and names of libraries are machine dependent

for WHICHIP in ctl test; do
  for PRECISION in 4 8 d; do

    case $PRECISION in
      4) SP_LIB=$SP_LIB4
         BACIO_LIB=$BACIO_LIB4
         W3NCO_LIB=$W3NCO_LIB4 ;;
      8) SP_LIB=$SP_LIB8
         BACIO_LIB=$BACIO_LIB8
         W3NCO_LIB=$W3NCO_LIB8 ;;
      d) SP_LIB=$SP_LIBd
         BACIO_LIB=$BACIO_LIB4
         W3NCO_LIB=$W3NCO_LIBd ;;
    esac

    ./configure --prefix=${PWD} --enable-promote=${PRECISION} \
      LIBS="${PWD}/lib/libip_${WHICHIP}_${PRECISION}.a ${SP_LIB} ${BACIO_LIB} ${W3NCO_LIB}"
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
