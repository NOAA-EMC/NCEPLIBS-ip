#!/bin/sh

set -x

icc -c -std=c99 test_gdswzd.c
icc test_gdswzd.o ../ip/v3.0.0/libip_v3.0.0_4.a /nwprod2/lib/sp/v2.0.2/libsp_v2.0.2_4.a /usrx/local/intel/composer_xe_2011_sp1.11.339/compiler/lib/intel64/libifcore.a -o test_gdswzd.exe
