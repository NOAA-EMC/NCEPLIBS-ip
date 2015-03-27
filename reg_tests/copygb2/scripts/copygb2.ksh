#!/bin/ksh

set -x

WORK_DIR=${WORK_DIR:-/stmpp1/$LOGNAME/regression}
WORK=${WORK_DIR}/copygb2
rm -fr $WORK
mkdir -p $WORK
cd $WORK

EXEC_DIR=/global/save/George.Gayno/iplib_branches/ip_grib2/reg_tests/copygb2/exec

# ncep grid 4 - regular lat/lon
grid[4]="0 6 255 -1 255 -1 255 -1 720 361 0 -1 90000000 0 56 -90000000 359500000 500000 500000 0"

# ncep grid 8 - mercator
grid[8]="10 6 255 -1 255 -1 255 -1 116 44 -48670000 3104000 56 22500000 61050000 0 64 0 318830000 318830000"

# ncep grid 127 - t254 gaussain
grid[127]="40 6 255 -1 255 -1 255 -1 768 384 0 -1 89642000 0 48 -89642000 359531000 469000 192 0"

# afwa grid 212 - nh polar stereographic
grid[212]="20 6 255 -1 255 -1 255 -1 512 512 -20826000 145000000 56 60000000 280000000 47625000 47625000 0 0"

# ncep grid 218 - lambert conformal
grid[218]="30 6 255 -1 255 -1 255 -1 614 428 12190000 226541000 56 25000000 265000000 12191000 12191000 0 64 25000000 25000000 -90000000 0"

for bytesize in 4 d 8 # test all three byte versions of library.
do
  COPYGB2=${EXEC_DIR}/copygb2_test_$bytesize
  for option in 0 1 2 3 "4 0 -1" 6
  do
    if [[ $option == "4 0 -1" ]];then # spectral option does not work with bitmaps.
      INPUT_FILE=/global/save/George.Gayno/iplib_branches/ip_grib2/reg_tests/copygb2/data/elevation.grb2
      opt=4
    else
      INPUT_FILE=/global/save/George.Gayno/iplib_branches/ip_grib2/reg_tests/copygb2/data/mxsnoalb.grb2
      opt=$option
    fi
    for gridnum in 4 8 127 212 218
    do
      $COPYGB2 -g "${grid[gridnum]}" -i"${option}" -x $INPUT_FILE grid${gridnum}.opt${opt}.${bytesize}byte.grb2
    done
  done
done

exit 0
