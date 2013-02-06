SHELL=  /bin/sh

MAKE_FILE=makefile

all:
	$(MAKE) -f $(MAKE_FILE) clobber
	$(MAKE) -f $(MAKE_FILE) 4_byte
	$(MAKE) -f $(MAKE_FILE) d_byte
	$(MAKE) -f $(MAKE_FILE) 8_byte

d_byte:
	export EXE=polateg_ctl_d.exe ; \
	if [ `uname -s` == "Linux" ];then  \
	  export FCOMP=ifort ; \
	  export FFLAGS="-check all -traceback -fpe0 -ftrapuv -assume byterecl -openmp -g -r8 -i4 -FR" ; \
	  export LIBIP="../../../lib/libip_ctl_d.a"  ; \
	  export LIBSP="/nwprod/lib/libsp_v2.0.1_d.a"  ; \
	  export LIBS="-L/nwprod/lib -lbacio_v2.0.1_4 -lw3nco_v2.0.4_d"  ; \
	elif [ `uname -s` == "AIX" ];then \
	  export FCOMP=xlf90_r ; \
	  export FFLAGS="-C -g -qrealsize=8 -qextchk -qsmp=omp" ; \
	  export LIBIP="../../../lib/libip_ctl_d.a"  ; \
	  export LIBS="-L/nwprod/lib -lbacio_4 -lw3_d -lsp_d -lessl"  ; \
	fi ; \
	$(MAKE) -f $(MAKE_FILE) build
8_byte:
	export EXE=polateg_ctl_8.exe ; \
	if [ `uname -s` == "Linux" ];then  \
	  export FCOMP=ifort ; \
	  export FFLAGS="-check all -traceback -fpe0 -ftrapuv -assume byterecl -openmp -g -r8 -i8 -FR" ; \
	  export LIBIP="../../../lib/libip_ctl_8.a"  ; \
	  export LIBSP="/nwprod/lib/libsp_v2.0.1_8.a"  ; \
	  export LIBS="-L/nwprod/lib -lbacio_v2.0.1_8 -lw3nco_v2.0.8_d"  ; \
	elif [ `uname -s` == "AIX" ];then \
	  export FCOMP=xlf90_r ; \
	  export FFLAGS="-C -g -qrealsize=8 -qintsize=8 -qextchk -qsmp=omp" ; \
	  export LIBIP="../../../lib/libip_ctl_8.a"  ; \
	  export LIBS="-L/nwprod/lib -lbacio_8 -lw3_8 -lsp_8 -lessl"  ; \
	fi ; \
	$(MAKE) -f $(MAKE_FILE) build
4_byte:
	export EXE=polateg_ctl_4.exe ; \
	if [ `uname -s` == "Linux" ];then  \
	  export FCOMP=ifort ; \
	  export FFLAGS="-check all -traceback -fpe0 -ftrapuv -assume byterecl -openmp -g -i4 -FR" ; \
	  export LIBIP="../../../lib/libip_ctl_4.a"  ; \
	  export LIBSP="/nwprod/lib/libsp_v2.0.1_4.a"  ; \
	  export LIBS="-L/nwprod/lib -lbacio_v2.0.1_4 -lw3nco_v2.0.4_d"  ; \
	elif [ `uname -s` == "AIX" ];then \
	  export FCOMP=xlf90_r ; \
	  export FFLAGS="-C -g -qextchk -qsmp=omp" ; \
	  export LIBIP="../../../lib/libip_ctl_4.a"  ; \
	  export LIBS="-L/nwprod/lib -lbacio_4 -lw3_4 -lsp_4 -lessl"  ; \
	fi ; \
	$(MAKE) -f $(MAKE_FILE) build

OBJS = driver.o \
       get_input_data.o  \
       interp.o

build : ${OBJS}
	$(FCOMP) $(FFLAGS) -o ${EXE} ${OBJS} $(LIBIP) $(LIBSP) $(LIBS)
	rm -f ${OBJS} *.mod

interp.o:
	${FCOMP} ${FFLAGS} -c ../../src/interp.f

get_input_data.o:
	${FCOMP} ${FFLAGS} -c ../../src/get_input_data.f

driver.o:  get_input_data.o interp.o
	${FCOMP} ${FFLAGS} -c ../../src/driver.f

clean :
	rm -f *.mod *.o 

clobber :
	rm -f *.mod *.o *.exe
