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
	  export FOPTS="-check all -traceback -fpe0 -ftrapuv -assume byterecl -g -r4 -i4" ; \
	  export LOPTS="-openmp" ; \
	  export LIBIP="../../../lib/libip_ctl_4.a" ; \
	  export LIBSP="/scratch1/portfolios/NCEPDEV/da/save/George.Gayno/sp_from_georgev/libsp_4.a" ; \
	  export LIBS="-L/scratch1/portfolios/NCEPDEV/da/save/Michael.Lueken/nwprod/lib -lw3lib-2.0_4 -lbacio_4" ; \
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
	  export LIBSP="/scratch1/portfolios/NCEPDEV/da/save/George.Gayno/sp_from_georgev/libsp_8.a" ; \
	  export LIBS="-L/scratch1/portfolios/NCEPDEV/da/save/Michael.Lueken/nwprod/lib -lw3lib-2.0_8 -lbacio_8" ; \
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
	  export LIBSP="/scratch1/portfolios/NCEPDEV/da/save/George.Gayno/sp_from_georgev/libsp_d.a" ; \
	  export LIBS="-L/scratch1/portfolios/NCEPDEV/da/save/Michael.Lueken/nwprod/lib -lw3lib-2.0_d -lbacio_4" ; \
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
