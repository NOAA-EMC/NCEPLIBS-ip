#!/bin/ksh --login

set -x

if [[ "$(hostname -f)" == g????.ncep.noaa.gov || \
      "$(hostname -f)" == t????.ncep.noaa.gov ]]; then  #WCOSS Phase 1/2
  module purge
  module load ics
  module load sp
fi

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
  icc -c -std=c99 -I../lib/incmod_${precision} test_gdswzd_${precision}.c
  icc test_gdswzd_${precision}.o ../lib/libip_${precision}.a ${SP_LIB} -lifcore -o test_gdswzd_${precision}.exe
  mv test_gdswzd_${precision}.exe ../exec
  rm -f *.o
done
