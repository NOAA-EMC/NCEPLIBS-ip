#!/bin/ksh

set -x

COPYGB=/nwprod/util/exec/copygb
infile=/nwprod/fix/global_shdmax.0.144x0.144.grb
outfile=./green.202.grb
kgds="255 202 955 835 -7491 -144134 136 955 835 126 90 64 0"

$COPYGB -g "$kgds" -i 2 -x $infile $outfile

exit 0
