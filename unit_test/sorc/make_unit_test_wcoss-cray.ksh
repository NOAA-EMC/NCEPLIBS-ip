#!/bin/ksh

#-----------------------------------------------------------------------------
# This script compiles the iplib unit test programs on the WCOSS-Cray
# system ONLY.  On all other machines, use "make_unit_test.ksh".
#
# PLEASE READ THE "README" FILE IN THIS DIRECTORY FOR DETAILS ON HOW
# TO RUN THIS SCRIPT.
#-----------------------------------------------------------------------------

usage()
{
 echo; echo "Usage: $0 setup-file" >&2
 echo; echo "Use this script on WCOSS-Cray only."
 echo; echo "Currently available setup files are:"
 echo
 for file in `ls ../../config-setup/crayftn* ` `ls ../../config-setup/ifort* ` \
             `ls ../../config-setup/gfortran* ` ; do
   echo "`basename ${file}`" >&2
 done
 echo
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
# The unit tests depend on the NCEP SP library.  The Intel and Cray 
# compiled versions of SP are found via modules.  Otherwise, these 
# environment variables must be set manually.
#-----------------------------------------------------------------------------

if [[ "$(hostname)" == slogin? || "$(hostname)" == llogin? ]]; then # WCOSS Cray ]]
  . /opt/modules/3.2.6.7/init/ksh
  module purge
  case $FC in
    ifort)
      FC="ftn"
      R8FLAG="-r8"
      I8FLAG="-i8"
      FCFLAGS="${FCFLAGS} -axCore-AVX2"  # add for haswell.
      module load PrgEnv-intel
      module load craype-sandybridge
      module load sp-intel/2.0.2 ;;
    crayftn)
      FC="ftn"
      R8FLAG="-s real64"
      I8FLAG="-s integer64"
      module load PrgEnv-cray
      module load craype-haswell
      module load sp-cray-haswell/2.0.2 ;;
    gfortran)
      FC="ftn"
      R8FLAG="-fdefault-real-8"
      I8FLAG="-fdefault-integer-8"
      module load PrgEnv-gnu
      module load craype-haswell ;;
  esac
else
  echo; echo "$0: ERROR - Script can only be used on WCOSS-Cray" >&2
  usage
  exit 13
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
# Make unit test executables for all three precision versions of IPLIB.
#-----------------------------------------------------------------------------

for PRECISION in 4 8 d; do  # single ("4"), double ("8") or mixed ("d") precison IPLIB

  case $PRECISION in
    4) SP_LIB=$SP_LIB4 
       FCFLAGS_ALL="${FCFLAGS}" ;;
    8) SP_LIB=$SP_LIB8 
       FCFLAGS_ALL="${FCFLAGS} ${R8FLAG} ${I8FLAG}" ;;
    d) SP_LIB=$SP_LIBd 
       FCFLAGS_ALL="${FCFLAGS} ${R8FLAG}" ;;
  esac

  ./configure --prefix=${root} --enable-promote=${PRECISION} \
    FC="${FC}" FCFLAGS="${FCFLAGS_ALL} -I${root}/lib/incmod_${PRECISION}" \
    LIBS="${root}/lib/libip_${PRECISION}.a ${SP_LIB}"
  if [ $? -ne 0 ]; then
    set +x
    echo "$0: Error configuring for ${PRECISION}-byte build." >&2
    exit 2
  fi

  $MAKE clean
  $MAKE all
  if [ $? -ne 0 ]; then
    set +x
    echo "$0: Error building ${PRECISION}-byte." >&2
    exit 3
  fi

  $MAKE install
  if [ $? -ne 0 ]; then
    set +x
    echo "$0: Error installing ${PRECISION}-byte." >&2
    exit 4
  fi

  mv config.log config_${PRECISION}.log

done  # library precision

echo DONE
