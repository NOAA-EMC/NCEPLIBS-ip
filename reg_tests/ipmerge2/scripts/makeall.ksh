#!/bin/ksh

set -x

cd ../exec/ctl
make all

cd ../test
make all

exit 0
