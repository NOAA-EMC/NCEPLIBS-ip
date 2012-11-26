SHELL=	/bin/sh
MAKE_FILE=makefile

all:
	make 4_byte
	make 8_byte
	make d_byte

4_byte:
	export CMD=copygb_ctl_4 ; \
	if [ `uname -s` == "Linux" ];then  \
	  export FCOMP="ifort" ; \
	  export FOPTS="-check all -traceback -fpe0 -ftrapuv -assume byterecl -g -i4" ; \
	  export LOPTS="-openmp" ; \
	  export LIBIP="../../../lib/libip_ctl_4.a" ; \
          export LIBSP="/nwprod/lib/libsp_v2.0.0_4.a"  ; \
          export LIBS="-L/nwprod/lib -lbacio_v2.0.1_4 -lw3nco_v2.0.3_4"  ; \
	elif [ `uname -s` == "AIX" ];then  \
	  export FCOMP="ncepxlf" ; \
	  export FOPTS="-O3 -qextchk  -qrealsize=4 -qnosave" ; \
	  export LOPTS="-qsmp=noauto  -lessl" ; \
	  export LIBIP="../../../lib/libip_ctl_4.a" ; \
	  export LIBSP="  " ; \
	  export LIBS="-L/nwprod/lib -lw3_4 -lbacio_4 -lsp_4" ; \
	fi ; \
	$(MAKE) -f $(MAKE_FILE) build

8_byte:
	export CMD=copygb_ctl_8 ; \
	if [ `uname -s` == "Linux" ];then  \
	  export FCOMP="ifort" ; \
	  export FOPTS="-check all -traceback -fpe0 -ftrapuv -assume byterecl -g -r8 -i8" ; \
	  export LOPTS="-openmp" ; \
	  export LIBIP="../../../lib/libip_ctl_8.a" ; \
          export LIBSP="/nwprod/lib/libsp_v2.0.0_8.a"  ; \
          export LIBS="-L/nwprod/lib -lbacio_v2.0.1_8 -lw3nco_v2.0.3_8"  ; \
	elif [ `uname -s` == "AIX" ];then  \
	  export FCOMP="ncepxlf" ; \
	  export FOPTS="-O3 -qextchk -qintsize=8 -qrealsize=8 -qnosave" ; \
	  export LOPTS="-qsmp=noauto  -lessl" ; \
	  export LIBIP="../../../lib/libip_ctl_8.a" ; \
	  export LIBSP="  " ; \
	  export LIBS="-L/nwprod/lib -lw3_8 -lbacio_8 -lsp_8" ; \
	fi ; \
	$(MAKE) -f $(MAKE_FILE) build

d_byte:
	export CMD=copygb_ctl_d ; \
	if [ `uname -s` == "Linux" ];then  \
	  export FCOMP="ifort" ; \
	  export FOPTS="-check all -traceback -fpe0 -ftrapuv -assume byterecl -g -r8 -i4" ; \
	  export LOPTS="-openmp" ; \
	  export LIBIP="../../../lib/libip_ctl_d.a" ; \
          export LIBSP="/nwprod/lib/libsp_v2.0.0_d.a"  ; \
          export LIBS="-L/nwprod/lib -lbacio_v2.0.1_4 -lw3nco_v2.0.3_d"  ; \
	elif [ `uname -s` == "AIX" ];then  \
	  export FCOMP="ncepxlf" ; \
	  export FOPTS="-O3 -qextchk -qintsize=4 -qrealsize=8 -qnosave" ; \
	  export LOPTS="-qsmp=noauto  -lessl" ; \
	  export LIBIP="../../../lib/libip_ctl_d.a" ; \
	  export LIBSP="  " ; \
	  export LIBS="-L/nwprod/lib -lw3_d -lbacio_4 -lsp_d" ; \
	fi ; \
	$(MAKE) -f $(MAKE_FILE) build

build:	copygb.o
	$(FCOMP) copygb.o $(LOPTS) $(LIBIP) $(LIBSP) $(LIBS) -o $(CMD)
	rm -f copygb.o

copygb.o:
	${FCOMP} -c $(FOPTS) ../../src/copygb.f
clean:	
	-rm -f copygb.o 
