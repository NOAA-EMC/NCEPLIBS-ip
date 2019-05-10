#!/bin/sh

#---------------------------------------------------------------
# Script to compile the regression test.
#
# Before invoking this script, the 'control' and 'test' libraries
# must be hosted in the ./reg_tests/lib directory with the
# following naming convention:
#
# libip_ctl_4.a  (4 byte integer/4 byte float, control)
# libip_ctl_8.a  (8 byte integer/8 byte float, control)
# libip_ctl_d.a  (4 byte integer/8 byte float, control)
# libip_test_4.a (4 byte integer/4 byte float, test)
# libip_test_8.a (8 byte integer/8 byte float, test)
# libip_test_d.a (4 byte integer/8 byte float, test)
#
# The 'control' and 'test' executables will be located in
# the ./exec/ctl and ./exec/test subdirectories.
# There will be one executable for each library version.
#---------------------------------------------------------------

set -x

cd ../exec/ctl 

make -f makefile clean
make -f makefile all

cd ../test

make -f makefile clean
make -f makefile all

exit 0
