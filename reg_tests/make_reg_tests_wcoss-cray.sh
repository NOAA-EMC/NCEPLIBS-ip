#!/bin/sh --login

#-----------------------------------------------------------------------------
# This script compiles all regression tests on the WCOSS-Cray machine only!
# DO NOT USE THIS SCRIPT ON OTHER MACHINES.
#
# $Id$
#
# To compile with the Intel compiler, uncomment the following line below:
#   . ./config-setup/ifort.setup
#
# To compile with the Cray compiler, uncomment the following line below:
#   . ./config-setup/crayftn.setup
#
# PLEASE READ THE "README" FILE IN THIS DIRECTORY FOR MORE DETAILS ON HOW
# TO RUN THIS SCRIPT.
#-----------------------------------------------------------------------------

#set -x

#-----------------------------------------------------------------------------
# These regression tests depend on the NCEP BACIO, SP, and W3NCO libraries.
# The path/name of these libraries are set thru modules.
#-----------------------------------------------------------------------------

if [[ "$(hostname)" == slogin? || "$(hostname)" == llogin? ]]; then # WCOSS Cray
  module purge
  module load modules/3.2.6.7
  module use /gpfs/hps/nco/ops/nwprod/lib/modulefiles
else
  echo
  echo "$0: Error. Script runs on WCOSS-Cray only. Abort." >&2
  exit 5
fi 

#-----------------------------------------------------------------------------
# Read in compiler, compiler flags and link flags.  
#-----------------------------------------------------------------------------

# uncomment to use Intel compiler
. ./config-setup/ifort.setup

# uncomment to use Cray compiler
#. ./config-setup/crayftn.setup

case $FC in
  ifort)
    module load PrgEnv-intel
    module load craype-haswell
    module load bacio-intel
    module load w3nco-intel
    module load sp-intel
    R8FLAG="-r8"
    I8FLAG="-i8"
    ;;
  crayftn)
    module load PrgEnv-cray
    module load craype-haswell
    module load sp-cray-haswell
    module load bacio-cray-haswell
    module load w3nco-cray-haswell
    R8FLAG="-s real64"
    I8FLAG="-s integer64"
    ;;
*)
  echo
  echo "$0: Error. Unrecognized compiler. Abort." >&2
  exit 8
  ;;
esac

module list

#-----------------------------------------------------------------------------
# Stop scripts if library environment variables are undefined.
#-----------------------------------------------------------------------------

BACIO_LIB4=${BACIO_LIB4:?}  # Single precision libraries
SP_LIB4=${SP_LIB4:?}
W3NCO_LIB4=${W3NCO_LIB4:?}

BACIO_LIB8=${BACIO_LIB8:?}  # Double precision libraries
SP_LIB8=${SP_LIB8:?}
W3NCO_LIB8=${W3NCO_LIB8:?}

SP_LIBd=${SP_LIBd:?}        # Mixed precision libraries
W3NCO_LIBd=${W3NCO_LIBd:?}

MAKE="gmake"

#-----------------------------------------------------------------------------
# Make regression test executables for all three precision versions of
# the 'control' and 'test IPLIB.
#-----------------------------------------------------------------------------

for WHICHIP in ctl test; do  # the 'control' or 'test' IPLIB
  for PRECISION in 4 8 d; do  # single ("4"), double ("8") or mixed ("d") precison IPLIB

    case $PRECISION in
      4) SP_LIB=$SP_LIB4
         BACIO_LIB=$BACIO_LIB4
         W3NCO_LIB=$W3NCO_LIB4
         FCFLAGS_ALL=${FCFLAGS} ;;
      8) SP_LIB=$SP_LIB8
         BACIO_LIB=$BACIO_LIB8
         W3NCO_LIB=$W3NCO_LIB8
         FCFLAGS_ALL="${FCFLAGS} ${R8FLAG} ${I8FLAG}" ;;
      d) SP_LIB=$SP_LIBd
         BACIO_LIB=$BACIO_LIB4
         W3NCO_LIB=$W3NCO_LIBd
         FCFLAGS_ALL="${FCFLAGS} ${R8FLAG}" ;;
    esac

    echo; echo
    echo "----------------------------------------------------------------------------------"
    echo "$0: Build ${PRECISION}-byte ${WHICHIP} version regression test suite."
    echo "----------------------------------------------------------------------------------"
    echo; echo

    ./configure --prefix=${PWD} --enable-promote=${PRECISION} \
      FC="ftn" FCFLAGS="${FCFLAGS_ALL} -I${PWD}/lib/incmod_${WHICHIP}_${PRECISION}" \
      LIBS="${PWD}/lib/libip_${WHICHIP}_${PRECISION}.a ${SP_LIB} ${BACIO_LIB} ${W3NCO_LIB}"
    if [ $? -ne 0 ]; then
      echo
      echo "$0: Error configuring for ${PRECISION}-byte ${WHICHIP} version build." >&2
      exit 2
    fi

    $MAKE clean
    $MAKE
    if [ $? -ne 0 ]; then
      echo
      echo "$0: Error building for ${PRECISION}-byte ${WHICHIP} version build." >&2
      exit 3
    fi

    $MAKE install suffix="_${WHICHIP}_${PRECISION}"
    if [ $? -ne 0 ]; then
      echo
      echo "$0: Error installing for ${PRECISION}-byte ${WHICHIP} version build." >&2
      exit 4
    fi

    mv config.log config_${WHICHIP}_${PRECISION}.log

  done  # library precision
done  # 'ctl' or 'test' IPLIB

echo; echo
echo "-------------------------------------------------------"
echo "$0: Done."
echo "-------------------------------------------------------"
