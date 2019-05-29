SHELL=  /bin/sh

MAKE_FILE=makefile

all:
	$(MAKE) -f $(MAKE_FILE) clobber
	$(MAKE) -f $(MAKE_FILE) 4_byte
	$(MAKE) -f $(MAKE_FILE) d_byte
	$(MAKE) -f $(MAKE_FILE) 8_byte

d_byte:
	export EXE=gausslat_ctl_d.exe ; \
	if [ `uname -s` == "Linux" ];then  \
	  export FCOMP=ifort ; \
	  export FFLAGS="-check all -traceback -fpe0 -ftrapuv -assume byterecl -g -r8 -i4 -FR" ; \
	  export LIBIP="../../../lib/libip_ctl_d.a"  ; \
	  export LIBSP="/u/Eugene.Mirvis/sp_merge/SP_MERGED/lib/libsp_d.a"  ; \
	  export LIBS="-L/nwprod/lib -lbacio_4 -lw3_d"  ; \
	elif [ `uname -s` == "AIX" ];then \
	  export FCOMP=xlf90_r ; \
	  export FFLAGS="-C -g -qrealsize=8 -qextchk" ; \
	  export LIBIP="../../../lib/libip_ctl_d.a"  ; \
	  export LIBS="-L/nwprod/lib -lbacio_4 -lw3_d -lsp_d -lessl"  ; \
	fi ; \
	$(MAKE) -f $(MAKE_FILE) build
8_byte:
	export EXE=gausslat_ctl_8.exe ; \
	if [ `uname -s` == "Linux" ];then  \
	  export FCOMP=ifort ; \
	  export FFLAGS="-check all -traceback -fpe0 -ftrapuv -assume byterecl -g -r8 -i8 -FR" ; \
	  export LIBIP="../../../lib/libip_ctl_8.a"  ; \
	  export LIBSP="/u/Eugene.Mirvis/sp_merge/SP_MERGED/lib/libsp_8.a"  ; \
	  export LIBS="-L/nwprod/lib -lbacio_8 -lw3_8"  ; \
	elif [ `uname -s` == "AIX" ];then \
	  export FCOMP=xlf90_r ; \
	  export FFLAGS="-C -g -qrealsize=8 -qintsize=8 -qextchk" ; \
	  export LIBIP="../../../lib/libip_ctl_8.a"  ; \
	  export LIBS="-L/nwprod/lib -lbacio_8 -lw3_8 -lsp_8 -lessl"  ; \
	fi ; \
	$(MAKE) -f $(MAKE_FILE) build
4_byte:
	export EXE=gausslat_ctl_4.exe ; \
	if [ `uname -s` == "Linux" ];then  \
	  export FCOMP=ifort ; \
	  export FFLAGS="-check all -traceback -fpe0 -ftrapuv -assume byterecl -g -i4 -FR" ; \
	  export LIBIP="../../../lib/libip_ctl_4.a"  ; \
	  export LIBSP="/u/Eugene.Mirvis/sp_merge/SP_MERGED/lib/libsp_4.a"  ; \
	  export LIBS="-L/nwprod/lib -lbacio_4 -lw3_4"  ; \
	elif [ `uname -s` == "AIX" ];then \
	  export FCOMP=xlf90_r ; \
	  export FFLAGS="-C -g -qextchk" ; \
	  export LIBIP="../../../lib/libip_ctl_4.a"  ; \
	  export LIBS="-L/nwprod/lib -lbacio_4 -lw3_4 -lsp_4 -lessl"  ; \
	fi ; \
	$(MAKE) -f $(MAKE_FILE) build

OBJS = gausslat_driver.o

build : ${OBJS}
	$(FCOMP) $(FFLAGS) -o ${EXE} ${OBJS} $(LIBIP) $(LIBSP) $(LIBS)
	rm -f ${OBJS}

gausslat_driver.o:
	${FCOMP} ${FFLAGS} -c ../../src/gausslat_driver.f

clean :
	rm -f *.o 

clobber :
	rm -f *.o *.exe
