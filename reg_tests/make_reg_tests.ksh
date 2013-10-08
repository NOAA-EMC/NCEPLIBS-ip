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
 
typeset -L1 MACHINE
MACHINE=$(hostname)

case $MACHINE in
# wcoss
  g|t) COMPILER="ifort" 
       COMPILER_FLAGS="-check all -traceback -fpe0 -ftrapuv -assume byterecl -g" 
       export COMPILER_FLAGS_XTRA="-openmp"
       export LD_FLAGS="-L${PWD}/lib -L/nwprod/lib" 
       export LD_FLAGS_XTRA="-openmp"
       SP="lsp_v2.0.1" 
       W3="lw3nco_v2.0.4"
       BACIO="lbacio_v2.0.1" ;;
# zeus
  f)   COMPILER="ifort" 
       COMPILER_FLAGS="-check all -traceback -fpe0 -ftrapuv -assume byterecl -g" 
       export COMPILER_FLAGS_XTRA="-openmp"
       export LD_FLAGS="-L${PWD}/lib -L/contrib/nceplibs/nwprod/lib" 
       export LD_FLAGS_XTRA="-openmp"
       SP="lsp_v2.0.1" 
       W3="lw3nco_v2.0.6"
       BACIO="lbacio_v2.0.1" ;;
# unknown machine
    *) set +x
       echo "$0: Error: Unknown Machine - Exiting" >&2
       exit 33 ;;
esac

MAKE="gmake"

# location and names of libraries are machine dependent

for WHICHIP in ctl test; do
  for PRECISION in 4 8 d; do

    case $PRECISION in
      d) PRECISION2=4 ;;
      *) PRECISION2=$PRECISION ;;
    esac

    ./configure --prefix=${PWD} --enable-promote=${PRECISION} FC=${COMPILER} FCFLAGS="${COMPILER_FLAGS}" \
      LDFLAGS="${LD_FLAGS}"  \
      LIBS="-lip_${WHICHIP}_${PRECISION} -${SP}_${PRECISION} -${BACIO}_${PRECISION2} -${W3}_${PRECISION}"
    ddd
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
