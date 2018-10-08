#!/bin/sh
###############################################################################
#
# $Id: make_ip_lib_wcoss-cray.sh 74394 2016-04-19 19:15:36Z george.gayno@noaa.gov $
#
# Script to iterate the configuration script over the set of precision
# versions of the library on the WCOSS-Cray.
#
# SCRIPT ONLY WORKS ON WCOSS-Cray.  Use make_ip_lib.sh on other machines.
#
# The build configuration setup (compiler, compiler switched, libraries, etc)
# is specified via files in the config-setup/ subdirectory that are sourced
# within this script.
#
# The installation directory is ${PWD}
#
###############################################################################

usage()
{
  echo
  echo " Usage: ${SCRIPT_NAME} [-g|-h] setup-file"
  echo
  echo "   Script to iterate the configuration script over the set of precision"
  echo "   versions of the library."
  echo
  echo "   ONLY RUN THIS SCRIPT ON WCOSS-CRAY!!"
  echo
  echo '   The installation directory is ${PWD}'
  echo
  echo " Options:"
  echo "   -g          Perform a Gnu-style install into include/ and lib/ directories"
  echo "               The default is an NCO-style install to reflect the structure"
  echo "               of /nwprod/lib"
  echo
  echo "   -h          Print this message and exit"
  echo
  echo " Arguments:"
  echo '   setup-file  File, in the "config-setup/" subdirectory, that contains'
  echo "               the build configuration setup (compiler, compiler switches,"
  echo "               libraries, etc) that are sourced within this script."
  echo
  echo "               Valid setup files for WCOSS-Cray are:"
  echo
  for file in `ls ./config-setup/crayftn* ./config-setup/ifort* ./config-setup/gfortran*`; do
    echo "     `basename ${file}`" >&2
  done
  echo
}

SCRIPT_NAME=$(basename $0)
SUCCESS=0
FAILURE=1
MAKE="gmake"
INSTALL_TYPE="nco"

if [[ "$(hostname)" != slogin? && "$(hostname)" != llogin? ]]; then # WCOSS Cray
  echo; echo "${SCRIPT_NAME}: ERROR - Script may only be run on WCOSS-Cray machine!"
  usage
  exit ${FAILURE}
fi

# Parse the command line options
while getopts :gh OPTVAL; do
  # Exit if option argument looks like another option
  case ${OPTARG} in
    -*) break;;
  esac
  # Parse the valid options
  case ${OPTVAL} in
    g) INSTALL_TYPE="gnu";;
    h)  usage
        exit ${SUCCESS};;
    :|\?) OPTVAL=${OPTARG}
          break;;
  esac
done
# ...Remove the options processed
shift $(expr ${OPTIND} - 1)
# ...Output invalidities based on OPTVAL
case ${OPTVAL} in
  # If OPTVAL contains nothing, then all options
  # have been successfully parsed and all that
  # remains are the arguments
  \?) if [ $# -lt 1 ]; then
        echo; echo "${SCRIPT_NAME}: ERROR - Missing build setup-file argument"
        usage
        exit ${FAILURE}
      fi;;
  # Invalid option
  ?) echo "${SCRIPT_NAME}: ERROR - Invalid option '-${OPTARG}'"
     usage
     exit ${FAILURE};;
esac

SETUP_FILE="./config-setup/$1"
if [ ! -f ${SETUP_FILE} ]; then
  echo
  echo "${SCRIPT_NAME}: ERROR - Cannot find specified setup file ${SETUP_FILE}" >&2
  usage
  exit ${FAILURE}
fi
. $SETUP_FILE

module purge
module load modules/3.2.6.7

case $FC in
  ifort)
    module load PrgEnv-intel
    module load craype-sandybridge     # Optimize for both sandybridge
    FCFLAGS="${FCFLAGS} -axCore-AVX2"  # and for haswell.
    R8FLAG="-r8"
    I8FLAG="-i8"
    FPPFLAGS="-fpp -save-temps"
    FPPFLAG4="-DLSIZE=4"
    FPPFLAGD="-DLSIZE=d"
    FPPFLAG8="-DLSIZE=8"
    COMP_NAME="intel"
    ;;
  gfortran)
    module load PrgEnv-gnu
    module load craype-haswell
    R8FLAG="-fdefault-real-8"
    I8FLAG="-fdefault-integer-8"
    FPPFLAGS="-cpp"
    FPPFLAG4="-DLSIZE=4"
    FPPFLAGD="-DLSIZE=d"
    FPPFLAG8="-DLSIZE=8"
    COMP_NAME="gnu"
    ;;
  crayftn)
    module load PrgEnv-cray
    module load craype-haswell
    R8FLAG="-s real64"
    I8FLAG="-s integer64"
    FPPFLAGS="-eZ"
    FPPFLAG4="-DLSIZE=4"
    FPPFLAGD="-DLSIZE=d"
    FPPFLAG8="-DLSIZE=8"
    COMP_NAME="cray"
    ;;
  *)
    echo "${SCRIPT_NAME}: ERROR - Unrecognized compiler ${FC}" >&2
    exit ${FAILURE} ;;
esac

echo
module list

# Build each precision version of library.

for PRECISION in 4 8 d; do  # single ("4"), double ("8") or mixed ("d") precison IPLIB

  case $PRECISION in
    4) FCFLAGS_ALL=${FCFLAGS} 
       FPPFLAGS_ALL="${FPPFLAGS} ${FPPFLAG4}";;
    8) FCFLAGS_ALL="${FCFLAGS} ${R8FLAG} ${I8FLAG}"
       FPPFLAGS_ALL="${FPPFLAGS} ${FPPFLAG8}";;
    d) FCFLAGS_ALL="${FCFLAGS} ${R8FLAG}"
       FPPFLAGS_ALL="${FPPFLAGS} ${FPPFLAGD}";;
  esac

  echo; echo; echo; echo
  echo "==============================================================="
  echo "==============================================================="
  echo "Configuring for precision ${PRECISION} build"
  echo "==============================================================="
  echo "==============================================================="
  echo

  ./configure --prefix=${PWD} --enable-promote=${PRECISION} --enable-wcoss_cray_dir=${COMP_NAME} \
    FC="ftn" FCFLAGS="${FCFLAGS_ALL} -craype-verbose" FPPFLAGS="${FPPFLAGS_ALL}"
  if [ $? -ne 0 ]; then
    echo "${SCRIPT_NAME}: ERROR configuring for precision ${PRECISION} version build" >&2
    exit ${FAILURE}
  fi

  echo; echo
  echo "==============================================================="
  echo "Starting precision ${PRECISION} build"
  echo "==============================================================="
  echo

  ${MAKE} clean
  ${MAKE}
  if [ $? -ne 0 ]; then
    echo "${SCRIPT_NAME}: ERROR building precision ${PRECISION} version" >&2
    exit ${FAILURE}
  fi

  # Install the current build...
  if [ "${INSTALL_TYPE}" = "nco" ]; then
    echo; echo
    echo "==============================================================="
    echo "Performing NCO-type install of precision ${PRECISION} build"
    echo "==============================================================="
    echo
    case $PRECISION in
      4) ${MAKE} nco_cray_uninstall 
    esac
    ${MAKE} nco_cray_install
    if [ $? -ne 0 ]; then
      echo "${SCRIPT_NAME}: ERROR in NCO-style installation of precision ${PRECISION} version" >&2
      exit ${FAILURE}
    fi
  else
    echo; echo
    echo "==============================================================="
    echo "Performing GNU-type install of precision ${PRECISION} build"
    echo "==============================================================="
    echo
    ${MAKE} uninstall
    ${MAKE} install
    if [ $? -ne 0 ]; then
      echo "${SCRIPT_NAME}: ERROR in Gnu-style installation of precision ${PRECISION} version" >&2
      exit ${FAILURE}
    fi
  fi

  # Clean up
  ${MAKE} distclean

done
