#!/bin/ksh

set -x
 
typeset -L1 MACHINE
MACHINE=$(hostname)

case $MACHINE in
# wcoss
  g|t) COMPILER="ifort" 
       FLAGS="-check all -traceback -fpe0 -ftrapuv -assume byterecl -g" 
       export OMP_FLAGS="-openmp"
       export LFLAGS="-L${PWD}/lib -L/nwprod/lib" 
       SP="lsp_v2.0.1" 
       W3="lw3nco_v2.0.4"
       BACIO="lbacio_v2.0.1" 
       EXTRA_LIB="" ;;
# zeus
  f)   COMPILER="ifort" 
       FLAGS="-check all -traceback -fpe0 -ftrapuv -assume byterecl -g -FR" 
       export OMP_FLAGS="-openmp"
       LFLAGS="-L${PWD}/lib -L/contrib/nceplibs/nwprod/lib" 
       SP="lsp_v2.0.1" 
       W3="lw3nco_v2.0.6"
       BACIO="lbacio_v2.0.1" 
       EXTRA_LIB="" ;;
# ccs
  c|s) COMPILER="xlf90_r" 
       FLAGS="-C -g -qextchk"
       export OMP_FLAGS="-qsmp=omp"
       LFLAGS="-L${PWD}/lib -L/nwprod/lib" 
       SP="lsp_v2.0.0" 
       W3="lw3_v2.2.3"
       BACIO="lbacio_v1.4.0"
       EXTRA_LIB="-lessl" ;;
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

    ./configure --prefix=${PWD} --enable-promote=${PRECISION} FC=${COMPILER} FCFLAGS="${FLAGS}" \
      LDFLAGS="${LFLAGS}"  \
      LIBS="-lip_${WHICHIP}_${PRECISION} -${SP}_${PRECISION} -${BACIO}_${PRECISION2} -${W3}_${PRECISION} ${EXTRA_LIB}"
    if [ $? -ne 0 ]; then
      set +x
      echo "$0: Error configuring for ${WHICHIP} precision ${PRECISION} version build" >&2
      exit 2
    fi

    $MAKE clean
    $MAKE
    if [ $? -ne 0 ]; then
      set +x
      echo "$0: Error building for ${WHICHIP} precision ${PRECISION} version build" >&2
      exit 3
    fi

    $MAKE install suffix="_${WHICHIP}_${PRECISION}"
    if [ $? -ne 0 ]; then
      set +x
      echo "$0: Error installing for ${WHICHIP} precision ${PRECISION} version build" >&2
     exit 4
    fi

    mv config.log config_${WHICHIP}_${PRECISION}.log

  done
done
