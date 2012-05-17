#!/bin/sh

set -x

cd ../exec/ctl 

make -f makefile clean
make -f makefile all

cd ../test

make -f makefile clean
make -f makefile all

exit 0
