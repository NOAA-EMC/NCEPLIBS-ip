#!/bin/sh

set -x

SETUP_FILE="./config-setup/$1"

. $SETUP_FILE

module purge
module load modules/3.2.6.7

case $FC in
  ifort)
    module load PrgEnv-intel
    R8FLAG="-r8"
    I8FLAG="-i8"
    ;;
  gfortran)
    module load PrgEnv-gnu
    R8FLAG="-fdefault-real-8"
    I8FLAG="-fdefault-integer-8"
    ;;
  crayftn)
    module load PrgEnv-cray
    R8FLAG="-s real64"
    I8FLAG="-s integer64"
    ;;
*)
  echo "unrecognized compiler"
  exit
  ;;
esac

module load craype-sandybridge
module list

# single precision version
./configure --prefix=${PWD} --enable-promote=4 FC="ftn"
make clean
make
make nco_install
make distclean

# mixed precision version
./configure --prefix=${PWD} --enable-promote=d FC="ftn" FCFLAGS="${R8FLAG} ${FCFLAGS}"
make clean
make
make nco_install
make distclean

# double precision version
./configure --prefix=${PWD} --enable-promote=8 FC="ftn" FCFLAGS="${I8FLAG} ${R8FLAG} ${FCFLAGS}"
make clean
make
make nco_install
make distclean

exit
