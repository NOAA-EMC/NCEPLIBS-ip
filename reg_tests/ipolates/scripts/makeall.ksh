#!/bin/ksh

set -x

cd ../exec/ctl
make clean
make all

cd ../test
make clean
make all

exit 0
