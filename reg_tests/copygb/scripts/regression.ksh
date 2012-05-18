#!/bin/ksh

# make sure i did not introduce errors into the ip and w3 libraries.

#set -x 

typeset -L4 machine
machine=$(dnsdomainname)

if [[ $machine == "zeus" ]];then    # zeus
  WORK_DIR=${WORK_DIR:-/scratch2/portfolios/NCEPDEV/stmp/$LOGNAME/regression}
else   # cirrus/stratus
  WORK_DIR=${WORK_DIR:-/stmp/$LOGNAME/regression}
fi

WORK=${WORK_DIR}/copygb
mkdir -p $WORK
cp ../exec/ctl/copygb_ctl_* $WORK
cp ../exec/test/copygb_test_* $WORK
cd $WORK

g2_test="./test.grb"
g2_ctl="./ctl.grb"
 
# file that will be interpolated
g1="/nwprod/fix/global_shdmax.0.144x0.144.grb"

for bytesize in 4 8 d
do
for option in 3  # interpolation option.  2 is nearest neighbor
                 # 3 is budget
do
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
  copygb_test_${bytesize} -g$grid -i${option} -x $g1 $g2_test
  copygb_ctl_${bytesize} -g$grid -i${option} -x $g1 $g2_ctl
  cmp test.grb ctl.grb
  status=$?
  if ((status != 0))
  then
    echo UNIX CMP FAILS FOR GRID $grid
    exit 8
  fi
  rm -f $g2_test $g2_ctl
done
done
done

exit 0
