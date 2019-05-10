#!/bin/ksh

#----------------------------------------------------------------------------
# Script to compile the entire suite of regression tests.
#
# Before invoking this script, the 'control' and 'test' libraries
# must be hosted in the ./reg_tests/lib directory with the following 
# naming convention:
#
# libip_ctl_4.a  (4 byte integer/4 byte float, control)
# libip_ctl_8.a  (8 byte integer/8 byte float, control)
# libip_ctl_d.a  (4 byte integer/8 byte float, control)
# libip_test_4.a (4 byte integer/4 byte float, test)
# libip_test_8.a (8 byte integer/8 byte float, test)
# libip_test_d.a (4 byte integer/8 byte float, test)
#
# Modify the links below for the libraries you wish to test.
#
# To invoke, type "Makeall.ksh"
#
# Six executables (one for each library) will be created for each
# regression test and be placed in the $testname/exec subdirectories.
#----------------------------------------------------------------------------

set -x

cd lib

#rm -f libip_ctl_4.a libip_ctl_8.a libip_ctl_d.a

#ln -fs /scratch1/portfolios/NCEPDEV/da/save/George.Gayno/iplib_g2/lib/libip_r19144_4.a  libip_ctl_4.a
#ln -fs /scratch1/portfolios/NCEPDEV/da/save/George.Gayno/iplib_g2/lib/libip_r19144_8.a  libip_ctl_8.a
#ln -fs /scratch1/portfolios/NCEPDEV/da/save/George.Gayno/iplib_g2/lib/libip_r19144_d.a  libip_ctl_d.a

#ln -fs /scratch1/portfolios/NCEPDEV/da/save/George.Gayno/iplib_g2/lib/libip_r19144_4.a  libip_test_4.a
#ln -fs /scratch1/portfolios/NCEPDEV/da/save/George.Gayno/iplib_g2/lib/libip_r19144_8.a  libip_test_8.a
#ln -fs /scratch1/portfolios/NCEPDEV/da/save/George.Gayno/iplib_g2/lib/libip_r19144_d.a  libip_test_d.a

cd ../copygb/scripts
makeall.ksh

cd ../../gausslat/scripts
makeall.ksh

cd ../../gcdist/scripts
makeall.ksh

cd ../../gdswiz_wzd/scripts
makeall.ksh

cd ../../ipmerge2/scripts
makeall.ksh

cd ../../ipolates/scripts
makeall.ksh

cd ../../ipolatev/scripts
makeall.ksh

cd ../../ipsector/scripts
makeall.ksh

cd ../../ipxetas/scripts
makeall.ksh

cd ../../ipxwafs/scripts
makeall.ksh

cd ../../makgds/scripts
makeall.ksh

exit 0
