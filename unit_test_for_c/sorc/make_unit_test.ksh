#!/bin/ksh --login

set -x

compiler=${compiler:-intel}

if [[ "$(hostname -f)" == g????.ncep.noaa.gov || \
      "$(hostname -f)" == t????.ncep.noaa.gov ]]; then  #WCOSS Phase 1/2
  case $compiler in
    intel)
      module purge
      module load ics
      module load sp ;;
    gnu)
      module purge
      CCOMP="gcc"
      CFLAGS="-std=c99"
      LIBS="-lgfortran" 
      SP_LIB4="/global/noscrub/George.Gayno/sp_v2.0.2/gfortran/libsp_v2.0.2_4.a"
      SP_LIB8="/global/noscrub/George.Gayno/sp_v2.0.2/gfortran/libsp_v2.0.2_8.a"
      SP_LIBd="/global/noscrub/George.Gayno/sp_v2.0.2/gfortran/libsp_v2.0.2_d.a" ;;
  esac
elif [[ "$(hostname)" == slogin? || "$(hostname)" == llogin? ]]; then # WCOSS Cray ]]
  . /opt/modules/3.2.6.7/init/ksh
  case $compiler in
    intel) module purge
           module load PrgEnv-intel
           module load craype-sandybridge
           module load sp-intel/2.0.2 
           CCOMP="cc"
           CFLAGS="-std=c99" 
           LIBS="-lifcore" ;; 
    cray)  module purge
           module load PrgEnv-cray
           module load craype-haswell
           module load sp-cray-haswell/2.0.2 
           CCOMP="cc"
           CFLAGS=" " 
           LIBS=" " ;;
  esac
elif [[ "$(hostname -f)" == tfe?? ]]; then # Theia
  case $compiler in
    intel)
      module purge
      module use -a /scratch3/NCEPDEV/nwprod/lib/modulefiles
      module load intel
      module load sp ;;
    gnu)
      module purge
      CCOMP="gcc"
      CFLAGS="-std=c99"
      LIBS="-lgfortran" 
      SP_LIB4="/scratch4/NCEPDEV/da/noscrub/George.Gayno/sp_v2.0.2/gfortran/libsp_v2.0.2_4.a"
      SP_LIB8="/scratch4/NCEPDEV/da/noscrub/George.Gayno/sp_v2.0.2/gfortran/libsp_v2.0.2_8.a"
      SP_LIBd="/scratch4/NCEPDEV/da/noscrub/George.Gayno/sp_v2.0.2/gfortran/libsp_v2.0.2_d.a" ;;
  esac
fi

CCOMP=${CCOMP:-icc}
CFLAGS=${CFLAGS:-"-std=c99"}
LIBS=${LIBS:-"-lifcore"}

SP_LIB4=${SP_LIB4:?}
SP_LIB8=${SP_LIB8:?}
SP_LIBd=${SP_LIBd:?}        # Mixed precision libraries

rm -f *.exe *.o
rm -f ../exec/*.exe

for precision in "4" "d" "8"
do
  case $precision in
    4) SP_LIB=$SP_LIB4 ;;
    8) SP_LIB=$SP_LIB8 ;;
    d) SP_LIB=$SP_LIBd ;;
  esac
  $CCOMP $CFLAGS -c -I../lib/incmod_${precision} test_gdswzd_${precision}.c
  $CCOMP test_gdswzd_${precision}.o ../lib/libip_${precision}.a ${SP_LIB} ${LIBS} -o test_gdswzd_${precision}.exe
  mv test_gdswzd_${precision}.exe ../exec
  rm -f *.o
done
