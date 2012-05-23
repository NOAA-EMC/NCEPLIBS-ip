#!/bin/ksh

#--------------------------------------------------------------------
# invoke copygb to interpolate a global lat/lon grid of vegetation
# greenness to numerous ncep standard grids.  use all interpolation
# options: 0-bilinear, 1-bicubic, 2-neighbor, 3-budget, 4-spectral,
# 6-neighbor-budget.
#
# the ncep standard grid numbers are defined in the official
# ncep grib documentation.
#
# If the output files from the test and control are not bit identical,
# the test fails.
#
# This is considered a 'supplimental' regression test and it not
# run as part of the full suite of tests from the "Runall" 
# driver script.
#--------------------------------------------------------------------

set -x 

WORK_DIR=${WORK_DIR:-/stmp/$LOGNAME/regression}

WORK=${WORK_DIR}/copygb
mkdir -p $WORK
cp ../exec/ctl/copygb_ctl_* $WORK
cp ../exec/test/copygb_test_* $WORK
cd $WORK

output_test="./test.grb"
output_ctl="./ctl.grb"
 
# file that will be interpolated
input="/nwprod/fix/global_shdmax.0.144x0.144.grb"

for bytesize in 4 8 d  # test all three byte versions of library.
do
  for option in 0 1 2 3 4 6  # interpolation option
  do     # ncep standard grid, see ncep grib 1 doc for details
    for grid in 1 2 3 4 5 6 8 10 11 12 13 14 15 16 17 18 \
             27 28 29 30 33 34 45 53 55 56 85 90 91 92 \
             93 94 95 96 97 98 99 100 101 103 104  \
             106 107 110 126 127 130 138 145 146 147 148 \
             150 151 160 161 163 170 171 172 173 174 175 176 \
             180 181 182 183 190 192 194 195 196 197 198 \
             201 202 203 204 205 206 207 208 209 210 211 212 \
             214 215 216 217 218 219 220 221 222 223 224 \
             225 226 227 228 229 230 231 232 233 234 235 236 237 238 239 \
             240 241 242 243 244 245 246 247 248 249 250 251 252 253 254
    do
      echo PROCESS GRID $grid
      copygb_test_${bytesize} -g$grid -i${option} -x $input $output_test
      copygb_ctl_${bytesize} -g$grid -i${option} -x $input $output_ctl
      cmp $output_test $output_ctl
      status=$?
      if ((status != 0))
      then
        echo UNIX CMP FAILS FOR GRID $grid
        exit 8
      fi
      rm -f $output_test $output_ctl
    done
  done
done

exit 0
