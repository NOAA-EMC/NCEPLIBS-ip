# *** manually set environments (for gnu compiler) of ip ***

 : ${USERMODE:=false}  # user mode (USERMODE) is closed by default
                       # set env var USERMODE to "true" to active it
 ${USERMODE} && {
    echo "Environment set by user"
# On theia/cray, user can load environment
#   module load gcc/6.2.0
# Or set environment on specific platform
    gcc_version=6.2.0
    gcc_topdir=/apps/gcc/$gcc_version
    export PATH=$gcc_topdir/bin:$PATH
    export LIBRARY_PATH=$gcc_topdir/lib64:$LIBRARY_PATH
    export LD_LIBRARY_PATH=$gcc_topdir/lib64:$LD_LIBRARY_PATH
    export INCLUDE=$gcc_topdir/include:$INCLUDE
    export MANPATH=$gcc_topdir/share/man:$MANPATH
 }

 ANCHORDIR=..
 export COMP=gnu
 export IP_VER=v2.0.0
 export IP_SRC=
 export IP_INC4=$ANCHORDIR/include/ip_${IP_VER}_4
 export IP_INC8=$ANCHORDIR/include/ip_${IP_VER}_8
 export IP_INCd=$ANCHORDIR/include/ip_${IP_VER}_d
 export IP_LIB4=$ANCHORDIR/libip_${IP_VER}_4.a
 export IP_LIB8=$ANCHORDIR/libip_${IP_VER}_8.a
 export IP_LIBd=$ANCHORDIR/libip_${IP_VER}_d.a

 export CC=gcc
 export FC=gfortran
 export CPP=cpp
 export OMPCC="$CC -fopenmp"
 export OMPFC="$FC -fopenmp"
 export MPICC=mpigcc
 export MPIFC=mpigfortran

 export DEBUG="-g -fbacktrace -O0"
 export CFLAGS="-g -O3 -fPIC"
 export FFLAGS="-g -fbacktrace -O3 -fconvert=little-endian -fPIC"
 export FREEFORM="-ffree-form"
 export FPPCPP="-cpp"
 export CPPFLAGS="-P -traditional-cpp"
 export MPICFLAGS="-g -O3 -fPIC"
 export MPIFFLAGS="-g -fbacktrace -O3 -fconvert=little-endian -fPIC"
 export MODPATH="-J"
 export I4R4=""
 export I4R8="-fdefault-real-8"
 export I8R8="-fdefault-integer-8 -fdefault-real-8"

 export CPPDEFS=""
 export CFLAGSDEFS="-DUNDERSCORE -DLINUX"
 export FFLAGSDEFS=""

 export USECC=""
 export USEFC="YES"
 export DEPS=""
