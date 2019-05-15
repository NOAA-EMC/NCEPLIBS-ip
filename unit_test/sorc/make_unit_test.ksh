#!/bin/ksh --login

#-----------------------------------------------------------------------------
# This script compiles the iplib unit test programs.
#
# PLEASE READ THE "README" FILE IN THIS DIRECTORY FOR DETAILS ON HOW
# TO RUN THIS SCRIPT.
#-----------------------------------------------------------------------------

usage()
{
 echo; echo "Usage: $0 setup-file" >&2
 echo; echo "Dont use this script on WCOSS-Cray."
 echo; echo "Currently available setup files are:"
 echo
 for file in `ls ../../config-setup/`; do
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
# The unit tests depend on the NCEP SP library.  Its path/name is set thru
# environment variables.  On Theia and WCOSS-Phase 1/2, the Intel
# version of SP is found thru modules.  On other machines, or when using
# a compiler other than Intel, the environment variables must be set
# manually.
#-----------------------------------------------------------------------------

if [[ "$(hostname -f)" == tfe?? ]]; then # Theia
  case $FC in
    ifort)
      module purge
      module use -a /scratch3/NCEPDEV/nwprod/lib/modulefiles
      module load intel/15.6.233
      module load sp ;;
  esac
elif [[ "$(hostname -f)" == g????.ncep.noaa.gov || \
        "$(hostname -f)" == t????.ncep.noaa.gov ]]; then  #WCOSS Phase 1/2
  case $FC in
    ifort)
      module purge
      module load ics/15.0.6
      module load sp ;;
  esac
elif [[ "$(hostname -f)" == v????.ncep.noaa.gov || \
        "$(hostname -f)" == m????.ncep.noaa.gov ]]; then  #WCOSS Dell
  case $FC in
    ifort)
      module purge
      module load EnvVars/1.0.2
      module load ips/18.0.1.163
      module load sp/2.0.2 ;;
  esac
elif [[ "$(hostname -f)" == slogin? || \
        "$(hostname -f)" == llogin? ]]; then  #WCOSS-Cray.
  echo; echo "$0: ERROR - Cant use this script on WCOSS-Cray." >&2
  usage
  exit 17
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
    4) SP_LIB=$SP_LIB4 ;;
    8) SP_LIB=$SP_LIB8 ;;
    d) SP_LIB=$SP_LIBd ;;
  esac

  ./configure --prefix=${root} --enable-promote=${PRECISION} \
    FC="${FC}" FCFLAGS="${FCFLAGS} -I${root}/lib/incmod_${PRECISION}" \
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
