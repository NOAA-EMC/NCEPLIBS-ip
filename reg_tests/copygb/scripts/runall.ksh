#!/bin/ksh

#--------------------------------------------------------------------
# Invoke copygb to interpolate a global lat/lon grid of vegetation
# greenness to numerous ncep standard grids.  Use all interpolation
# options: 0-bilinear, 1-bicubic, 2-neighbor, 3-budget, 4-spectral,
# 6-neighbor-budget.
#
# The ncep standard grid numbers are defined in the official
# ncep grib documentation.  Not all standard grids are tested. 
# For example, there is an error with how W3 routine w3fi71 
# specifies the GDS for NCEP grid 174.  This causes copygb to
# fail on Zeus. Also grids 90, 91 and 92 are very high-res, which
# will cause copygb to run very slow, especially for the spectral
# interpolation.  Running the spectral test just for grid 90
# took 44 minutes on Zeus.  
#
# If the output files from the test and control are not bit identical,
# the regression test fails.
#
# This is considered a 'supplemental' regression test and it not
# run as part of the full suite of tests from the main "Runall" 
# driver script.
#--------------------------------------------------------------------

#set -x 

echo
echo "BEGIN COPYGB REGRESSION TEST"
echo

REG_DIR=${REG_DIR:-../..}
EXEC_DIR=$REG_DIR/copygb/exec

WORK_DIR=${WORK_DIR:-/stmp/$LOGNAME/regression}
WORK=${WORK_DIR}/copygb
rm -fr $WORK
mkdir -p $WORK
WORK_CTL=${WORK}/ctl
mkdir -p $WORK_CTL
WORK_TEST=${WORK}/test
mkdir -p $WORK_TEST

cp $EXEC_DIR/ctl/copygb_ctl_* $WORK_CTL
cp $EXEC_DIR/test/copygb_test_* $WORK_TEST

output_ctl=$WORK_CTL/ctl.grb
output_test=$WORK_TEST/test.grb

reg_test_failed=0

for bytesize in 4 8 d  # test all three byte versions of library.
do

  echo TEST $bytesize BYTE VERSION OF LIBRARY.

  for option in "0" "1" "2" "3" "4 0 -1" "6"  # interpolation option
  do  

# The spectral interpolation option can run very slow. To speed things
# up, use a lower-res input data.

    if [[ $option == "4 0 -1" ]];then
      input="/nwprod/fix/global_snoalb.1x1.grb"
    else
      input="/nwprod/fix/global_shdmax.0.144x0.144.grb"
    fi

   # ncep standard grid, see ncep grib 1 doc for details
    for grid in 1 2 3 4 5 6 8 10 11 12 13 14 15 16 17 18 \
             27 28 29 30 33 34 45 53 55 56 85 \
             93 94 95 96 97 98 99 100 101 103 104  \
             106 107 110 126 127 130 138 145 146 147 148 \
             150 151 160 161 163 170 171 172 175 176 \
             180 181 182 183 190 192 194 195 196 197 198 \
             201 202 203 204 205 206 207 208 209 210 211 212 \
             214 215 216 217 218 219 220 221 222 223 224 \
             225 227 228 229 230 231 232 233 234 235 236 237 238 239 \
             240 241 242 243 244 245 246 247 248 249 250 251 252 253 254
    do

      echo PROCESS GRID ${grid} USING INTERPOLATION OPTION ${option}

      ctl_failed=0
      test_failed=0

      cd $WORK_CTL
      copygb_ctl_${bytesize} -g$grid -i"${option}" -x $input $output_ctl
      status=$?
      if ((status != 0)); then
        echo "** PROBLEM WITH CONTROL RUN."
        reg_test_failed=1
        ctl_failed=1
      fi

      if [[ ! -s $output_ctl ]];then
        echo "** CONTROL RUN CREATED ZERO BYTE FILE"
        reg_test_failed=1
        ctl_failed=1
      fi

      cd $WORK_TEST
      copygb_test_${bytesize} -g$grid -i"${option}" -x $input $output_test
      status=$?
      if ((status != 0)); then
        echo "** PROBLEM WITH TEST RUN."
        reg_test_failed=1
        test_failed=1
      fi

      if [[ ! -s $output_test ]];then
        echo "** TEST RUN CREATED ZERO BYTE FILE"
        reg_test_failed=1
        test_failed=1
      fi

      if ((ctl_failed == 0 && test_failed == 0));then
        cmp $output_test $output_ctl
        status=$?
        if ((status != 0)); then
          echo "** GRIB FILES NOT BIT IDENTICAL. REGRESSION TEST FAILS."
          reg_test_failed=1
        fi
      fi

      rm -f $output_test $output_ctl

    done
  done
done

if ((reg_test_failed == 0)); then
  echo
  echo "<<< COPYGB REGRESSION TEST PASSED. >>>"
  echo
else
  echo
  echo "<<< COPYGB REGRESSION TEST FAILED. >>>"
  echo
fi

exit 0
