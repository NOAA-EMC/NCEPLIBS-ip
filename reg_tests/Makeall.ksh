#!/bin/ksh

set -x

cd lib

rm -f libip_ctl_4.a libip_ctl_8.a libip_ctl_d.a

ln -fs /scratch1/portfolios/NCEPDEV/da/save/George.Gayno/iplib_g2/lib/libip_r19144_4.a  libip_ctl_4.a
ln -fs /scratch1/portfolios/NCEPDEV/da/save/George.Gayno/iplib_g2/lib/libip_r19144_8.a  libip_ctl_8.a
ln -fs /scratch1/portfolios/NCEPDEV/da/save/George.Gayno/iplib_g2/lib/libip_r19144_d.a  libip_ctl_d.a

ln -fs /scratch1/portfolios/NCEPDEV/da/save/George.Gayno/iplib_g2/lib/libip_r19144_4.a  libip_test_4.a
ln -fs /scratch1/portfolios/NCEPDEV/da/save/George.Gayno/iplib_g2/lib/libip_r19144_8.a  libip_test_8.a
ln -fs /scratch1/portfolios/NCEPDEV/da/save/George.Gayno/iplib_g2/lib/libip_r19144_d.a  libip_test_d.a
exit 0

cd ../reg_tests/copygb/scripts
makeall.ksh

cd ../gausslat/scripts
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

#cd ../polateg
#makeall.ksh
