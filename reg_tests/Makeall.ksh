#!/bin/ksh

set -x
 
# Check that a build setup file was specified.
if [ $# -ne 1 ]; then
  echo "$0: Must specify compiler." >&2
  exit 1
fi

COMPILER=$1

for whichip in ctl test
do

for PRECISION in 4 8 d
do

  case $PRECISION in
    d) PRECISION2=4 ;;
    *) PRECISION2=$PRECISION ;;
  esac

  ./configure --prefix=${PWD} --enable-promote=${PRECISION} FC=$COMPILER FCFLAGS="-openmp -FR" \
    LDFLAGS="-L${PWD}/lib -L/nwprod/lib" \
    LIBS="-lip_${whichip}_${PRECISION} -lsp_${PRECISION} -lbacio_${PRECISION2} -lw3nco_${PRECISION}"
  if [ $? -ne 0 ]; then
    echo "$0: Error configuring for ${whichip} precision ${PRECISION} version build" >&2
    exit 2
  fi
  gmake clean
  if [ $? -ne 0 ]; then
    echo "error"
    exit 7
  fi
  gmake
  if [ $? -ne 0 ]; then
    echo "error"
    exit 3
  fi
  gmake install suffix="_${whichip}_${PRECISION}"
  if [ $? -ne 0 ]; then
    echo "error"
    exit 4
  fi

  mv config.log config_${whichip}_${PRECISION}.log

done
done
