# *** manually set environments (for gnu compiler) of ip ***

# !!! module environment (*THEIA*) !!!
 module load gcc/6.2.0

 ANCHORDIR=..
 export COMP=gnu
 export IP_VER=v3.0.1
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

 export DEBUG="-g -O0"
 export CFLAGS="-O3 -DUNDERSCORE -DLINUX -fPIC"
 export FFLAGS="-O3 -fconvert=little-endian -fPIC"
 export CPPFLAGS="-P -traditional-cpp"
 export MPICFLAGS="-O3 -DUNDERSCORE -DLINUX -fPIC"
 export MPIFFLAGS="-O3 -fPIC"
 export MODPATH="-J"
 export I4R4=""
 export I4R8="-fdefault-real-8"
 export I8R8="-fdefault-integer-8 -fdefault-real-8"

 export CPPDEFS=""
 export CFLAGSDEFS=""
 export FFLAGSDEFS=""

 export USECC=""
 export USEFC="YES"
 export DEPS=""
