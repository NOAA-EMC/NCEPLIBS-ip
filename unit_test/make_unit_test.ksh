#!/bin/ksh --login

#-----------------------------------------------------------------------------
# This script compiles the iplib unit test.  The test uses the
# NCEP "copygb" and "diffgb" programs.
# 
# PLEASE READ THE "README" FILE IN THIS DIRECTORY FOR DETAILS ON HOW
# TO RUN THIS SCRIPT.
#-----------------------------------------------------------------------------

set -x

#-----------------------------------------------------------------------------
# The unit tests depend on the NCEP BACIO, SP, and W3NCO libraries.
# The path/name of these libraries are set thru environment variables.
# On Theia and WCOSS, these are set via modules.  On other machines,
# they must be set manually.
#-----------------------------------------------------------------------------

if [[ "$(hostname -f)" == tfe?? ]]; then # Theia
  module purge
  module use -a /scratch3/NCEPDEV/nwprod/lib/modulefiles
  module load intel
  module load bacio
  module load sp
  module load w3nco
elif [[ "$(hostname -d)" == "ncep.noaa.gov" ]]; then  #WCOSS
  module purge
  module load ics
  module load bacio
  module load sp
  module load w3nco
fi 

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

#-----------------------------------------------------------------------------
# Set some parameters.
#-----------------------------------------------------------------------------

MAKE="gmake"

root=${PWD}

#-----------------------------------------------------------------------------
# Make the "diffgb" component.
#-----------------------------------------------------------------------------

cd util

. ./config-setup/ifort.setup

./configure --prefix=${root} FCFLAGS="${FCFLAGS}" \
  LIBS="${SP_LIBd} ${BACIO_LIB4} ${W3NCO_LIBd}"
if [ $? -ne 0 ]; then
  set +x
  echo "$0: Error configuring for diffgb build." >&2
  exit 12
fi

$MAKE clean
$MAKE
if [ $? -ne 0 ]; then
  set +x
  echo "$0: Error building diffgb program." >&2
  exit 13
fi

$MAKE install
if [ $? -ne 0 ]; then
  set +x
  echo "$0: Error installing diffgb program." >&2
  exit 14
fi

#-----------------------------------------------------------------------------
# Make copygb executables for all three precision versions of IPLIB.
#-----------------------------------------------------------------------------

cd ../sorc

. ./config-setup/ifort.setup

for PRECISION in 4 8 d; do  # single ("4"), double ("8") or mixed ("d") precison IPLIB

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

  ./configure --prefix=${root} --enable-promote=${PRECISION} \
    FCFLAGS="${FCFLAGS} -I../lib/incmod_${PRECISION}" \
    LIBS="../lib/libip_${PRECISION}.a ${SP_LIB} ${BACIO_LIB} ${W3NCO_LIB}"
  if [ $? -ne 0 ]; then
    set +x
    echo "$0: Error configuring for ${PRECISION}-byte copygb build." >&2
    exit 2
  fi

  $MAKE clean
  $MAKE
  if [ $? -ne 0 ]; then
    set +x
    echo "$0: Error building ${PRECISION}-byte copygb." >&2
    exit 3
  fi

  $MAKE install
  if [ $? -ne 0 ]; then
    set +x
    echo "$0: Error installing ${PRECISION}-byte copygb." >&2
    exit 4
  fi

  mv config.log config_${PRECISION}.log

done  # library precision

echo DONE
