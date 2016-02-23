#!/bin/ksh --login

#-----------------------------------------------------------------------------
# This script compiles the iplib unit test.  The test uses the
# NCEP "copygb" and "diffgb" programs.
# 
# PLEASE READ THE "README" FILE IN THIS DIRECTORY FOR DETAILS ON HOW
# TO RUN THIS SCRIPT.
#-----------------------------------------------------------------------------

usage()
{
 echo; echo "Usage: $0 setup-file" >&2
 echo; echo "Currently available setup files are:"
 echo
 for file in `ls ../../config-setup/`; do
   echo "`basename ${file}`" >&2
 done
}

#set -x

#-----------------------------------------------------------------------------
# Script requires one argument - the name of the build setup file.
#-----------------------------------------------------------------------------

if [ $# -lt 1 ]; then
  echo; echo "$0: ERROR - Missing build setup file argument" >&2
  usage
  exit 19
fi

#-----------------------------------------------------------------------------
# Source the build setup
#-----------------------------------------------------------------------------

SETUP_FILE="../../config-setup/$1"
if [ ! -f ${SETUP_FILE} ]; then
  echo; echo "$0: ERROR - Cannot find specified setup file ${SETUP_FILE}" >&2
  usage
  exit 9
fi
. ${SETUP_FILE}

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
  module load sp
elif [[ "$(hostname -f)" == g????.ncep.noaa.gov || \
        "$(hostname -f)" == t????.ncep.noaa.gov ]]; then  #WCOSS
  case $FC in
    ifort)
      module purge
      module load ics
      module load sp ;;
  esac
fi 

#-----------------------------------------------------------------------------
# Stop scripts if library environment variables are undefined.
#-----------------------------------------------------------------------------

SP_LIB4=${SP_LIB4:?}
SP_LIB8=${SP_LIB8:?}
SP_LIBd=${SP_LIBd:?}        # Mixed precision libraries

#-----------------------------------------------------------------------------
# Set some parameters.
#-----------------------------------------------------------------------------

MAKE="gmake"

root="$PWD/.."

#-----------------------------------------------------------------------------
# Make copygb executables for all three precision versions of IPLIB.
#-----------------------------------------------------------------------------

for PRECISION in 4 8 d; do  # single ("4"), double ("8") or mixed ("d") precison IPLIB

  case $PRECISION in
    4) SP_LIB=$SP_LIB4 ;;
    8) SP_LIB=$SP_LIB8 ;;
    d) SP_LIB=$SP_LIBd ;;
  esac

  ./configure --prefix=${root} --enable-promote=${PRECISION} \
    FC="${FC}" FCFLAGS="${FCFLAGS} -I${root}/lib/incmod_${PRECISION}" \
    LIBS="${root}/lib/libip_${PRECISION}.a ${SP_LIB} ${BACIO_LIB} ${W3NCO_LIB}"
  if [ $? -ne 0 ]; then
    set +x
    echo "$0: Error configuring for ${PRECISION}-byte copygb build." >&2
    exit 2
  fi

  $MAKE clean
  $MAKE all
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
