SHELL=  /bin/sh

MAKE_FILE=makefile

all:
	$(MAKE) -f $(MAKE_FILE) clobber
	$(MAKE) -f $(MAKE_FILE) 4_byte
	$(MAKE) -f $(MAKE_FILE) d_byte
	$(MAKE) -f $(MAKE_FILE) 8_byte

d_byte:
	export EXE=ipolatev_test_d.exe ; \
	if [ `uname -s` == "Linux" ];then  \
	  export FCOMP=ifort ; \
	  export FFLAGS="-check all -traceback -fpe0 -ftrapuv -assume byterecl -g -r8 -i4 -FR -openmp" ; \
	  export LIBIP="../../../lib/libip_test_d.a"  ; \
	  export LIBSP="/scratch1/portfolios/NCEPDEV/da/save/George.Gayno/sp_from_georgev/libsp_d.a"  ; \
	  export LIBS="-L/scratch1/portfolios/NCEPDEV/da/save/Michael.Lueken/nwprod/lib -lbacio_4 -lw3lib-2.0_d"  ; \
	elif [ `uname -s` == "AIX" ];then \
	  export FCOMP=xlf90_r ; \
	  export FFLAGS="-C -g -qrealsize=8 -qextchk -qsmp=omp" ; \
	  export LIBIP="../../../lib/libip_test_d.a"  ; \
	  export LIBS="-L/nwprod/lib -lbacio_4 -lw3_d -lsp_d -lessl"  ; \
	fi ; \
	$(MAKE) -f $(MAKE_FILE) build
8_byte:
	export EXE=ipolatev_test_8.exe ; \
	if [ `uname -s` == "Linux" ];then  \
	  export FCOMP=ifort ; \
	  export FFLAGS="-check all -traceback -fpe0 -ftrapuv -assume byterecl -g -r8 -i8 -FR -openmp" ; \
	  export LIBIP="../../../lib/libip_test_8.a"  ; \
	  export LIBSP="/scratch1/portfolios/NCEPDEV/da/save/George.Gayno/sp_from_georgev/libsp_8.a"  ; \
	  export LIBS="-L/scratch1/portfolios/NCEPDEV/da/save/Michael.Lueken/nwprod/lib -lbacio_8 -lw3lib-2.0_8"  ; \
	elif [ `uname -s` == "AIX" ];then \
	  export FCOMP=xlf90_r ; \
	  export FFLAGS="-C -g -qrealsize=8 -qintsize=8 -qextchk -qsmp=omp" ; \
	  export LIBIP="../../../lib/libip_test_8.a"  ; \
	  export LIBS="-L/nwprod/lib -lbacio_8 -lw3_8 -lsp_8 -lessl"  ; \
	fi ; \
	$(MAKE) -f $(MAKE_FILE) build
4_byte:
	export EXE=ipolatev_test_4.exe ; \
	if [ `uname -s` == "Linux" ];then  \
	  export FCOMP=ifort ; \
	  export FFLAGS="-check all -traceback -fpe0 -ftrapuv -assume byterecl -g -i4 -FR -openmp" ; \
	  export LIBIP="../../../lib/libip_test_4.a"  ; \
	  export LIBSP="/scratch1/portfolios/NCEPDEV/da/save/George.Gayno/sp_from_georgev/libsp_4.a"  ; \
	  export LIBS="-L/scratch1/portfolios/NCEPDEV/da/save/Michael.Lueken/nwprod/lib -lbacio_4 -lw3lib-2.0_4"  ; \
	elif [ `uname -s` == "AIX" ];then \
	  export FCOMP=xlf90_r ; \
	  export FFLAGS="-C -g -qextchk -qsmp=omp" ; \
	  export LIBIP="../../../lib/libip_test_4.a"  ; \
	  export LIBS="-L/nwprod/lib -lbacio_4 -lw3_4 -lsp_4 -lessl"  ; \
	fi ; \
	$(MAKE) -f $(MAKE_FILE) build

OBJS = ipolatev_driver.o \
       get_input_data.o  \
       interp.o

build : ${OBJS}
	$(FCOMP) $(FFLAGS) -o ${EXE} ${OBJS} $(LIBIP) $(LIBSP) $(LIBS)
	rm -f ${OBJS} *.mod

interp.o:
	${FCOMP} ${FFLAGS} -c ../../src/interp.f

get_input_data.o:
	${FCOMP} ${FFLAGS} -c ../../src/get_input_data.f

ipolatev_driver.o:  get_input_data.o interp.o
	${FCOMP} ${FFLAGS} -c ../../src/ipolatev_driver.f

clean :
	rm -f *.mod *.o 

clobber :
	rm -f *.mod *.o *.exe
