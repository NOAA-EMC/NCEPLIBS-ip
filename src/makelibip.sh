#!/bin/sh
###############################################################
#
#   AUTHOR:    Vuong - W/NP11
#
#   DATE:      28/11/2000
#
#   PURPOSE:   This script uses the make utility to update the libip 
#              archive libraries.
#              It first reads a list of source files in the library and
#              then generates a makefile used to update the archive
#              libraries.  The make command is then executed for each
#              archive library, where the archive library name and 
#              compilation flags are passed to the makefile through 
#              environment variables.
#
#   REMARKS:   Only source files that have been modified since the last
#              library update are recompiled and replaced in the object
#              archive libraries.  The make utility determines this
#              from the file modification times.
#
#              New source files are also compiled and added to the object 
#              archive libraries.
#
###############################################################

#
#     Generate a list of object files that corresponds to the
#     list of Fortran ( .f90 ) files in the current directory
#
for i in `ls *.f`
do
  obj=`basename $i .f`
  OBJS="$OBJS ${obj}.o"
done
#
#     Remove make file, if it exists.  May need a new make file
#     with an updated object file list.
#
if [ -f make.libip ] 
then
  rm -f make.libip
fi
#
#     Generate a new make file ( make.libip), with the updated object list,
#     from this HERE file.
#
cat > make.libip << EOF
SHELL=/bin/sh

.SUFFIXES:
.SUFFIXES: .f90 .f .a

\$(LIB):	\$(LIB)( ${OBJS} )

.f90.a:
	\$(FCOMP) -c \$(FFLAGS) \$<
	ar -ruv \$(AFLAGS) \$@ \$*.o
	rm -f \$*.o

.f.a:
	\$(FCOMP) -c \$(FFLAGS) \$<
	ar -ruv \$(AFLAGS) \$@ \$*.o
	rm -f \$*.o

EOF
#
#     Update 4-byte version of libip_4.a
#
export LIB="../lib/libip_4.a"
if [ `uname -s` == "Linux" ];then
  export FCOMP="ifort"
  export FFLAGS="-check all -traceback -fpe0 -ftrapuv -g -r4 -i4 -openmp"
  export AFLAGS=" "
elif [ `uname -s` == "AIX" ];then
  export FCOMP="xlf90_r"
  export FFLAGS=" -O3 -qsmp=noauto -qnosave -qfixed"
  export AFLAGS=" -X64"
fi
make -f make.libip
#
#     Update 8-byte version of libip_8.a
#
export LIB="../lib/libip_8.a"
if [ `uname -s` == "Linux" ];then
  export FCOMP="ifort"
  export FFLAGS="-check all -traceback -fpe0 -ftrapuv -g -r8 -i8 -openmp"
  export AFLAGS=" "
elif [ `uname -s` == "AIX" ];then
  export FCOMP="xlf90_r"
  export FFLAGS=" -O3 -qsmp=noauto -qnosave -qintsize=8 -qrealsize=8 -qfixed"
  export AFLAGS=" -X64"
fi
make -f make.libip
#
#     Update Double Precision (Size of Real 8-byte and default Integer) version
#     of libip_d.a
#
export LIB="../lib/libip_d.a"
if [ `uname -s` == "Linux" ];then
  export FCOMP="ifort"
  export FFLAGS="-check all -traceback -fpe0 -ftrapuv -g -r8 -i4 -openmp"
  export AFLAGS=" "
elif [ `uname -s` == "AIX" ];then
  export FCOMP="xlf90_r"
  export FFLAGS=" -O3 -qsmp=noauto -qnosave -qintsize=4 -qrealsize=8 -qfixed"
  export AFLAGS=" -X64"
fi
make -f make.libip
rm -f make.libip
