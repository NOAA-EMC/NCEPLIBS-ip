#!/bin/sh

set -x

rm -f *.exe *.o

icc -c -std=c99 test_gdswzd_4.c
icc test_gdswzd_4.o ../ip/v3.0.0/libip_v3.0.0_4.a /nwprod2/lib/sp/v2.0.2/libsp_v2.0.2_4.a /usrx/local/intel/composer_xe_2011_sp1.11.339/compiler/lib/intel64/libifcore.a -o test_gdswzd_4.exe

icc -c -std=c99 test_gdswzd_d.c
icc test_gdswzd_d.o ../ip/v3.0.0/libip_v3.0.0_d.a /nwprod2/lib/sp/v2.0.2/libsp_v2.0.2_d.a /usrx/local/intel/composer_xe_2011_sp1.11.339/compiler/lib/intel64/libifcore.a -o test_gdswzd_d.exe

icc -c -std=c99 test_gdswzd_8.c
icc test_gdswzd_8.o ../ip/v3.0.0/libip_v3.0.0_8.a /nwprod2/lib/sp/v2.0.2/libsp_v2.0.2_8.a /usrx/local/intel/composer_xe_2011_sp1.11.339/compiler/lib/intel64/libifcore.a -o test_gdswzd_8.exe

rm -f *.o
