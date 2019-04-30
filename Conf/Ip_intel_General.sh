# *** manually set environments (for intel compiler) of ip ***

# !!! module environment (*THEIA*) !!!
 module load intel/18.1.163
#module load ics/17.0.3

 ANCHORDIR=..
 export COMP=ips
 export IP_VER=v3.0.1
 export IP_SRC=
 export IP_INC4=$ANCHORDIR/include/ip_${IP_VER}_4
 export IP_INC8=$ANCHORDIR/include/ip_${IP_VER}_8
 export IP_INCd=$ANCHORDIR/include/ip_${IP_VER}_d
 export IP_LIB4=$ANCHORDIR/libip_${IP_VER}_4.a
 export IP_LIB8=$ANCHORDIR/libip_${IP_VER}_8.a
 export IP_LIBd=$ANCHORDIR/libip_${IP_VER}_d.a

 export CC=icc
 export FC=ifort
 export CPP=cpp
 export OMPCC="$CC -qopenmp"
 export OMPFC="$FC -qopenmp"
 export MPICC=mpiicc
 export MPIFC=mpiifort

 export DEBUG="-g -O0"
 export CFLAGS="-O3 -DUNDERSCORE -DLINUX -fPIC"
 export FFLAGS="-O3 -fp-model strict -ip -convert little_endian -assume byterecl -fPIC"
 export CPPFLAGS="-P -traditional-cpp"
 export MPICFLAGS="-O3 -DUNDERSCORE -DLINUX -fPIC"
 export MPIFFLAGS="-O3 -xHOST -traceback -fPIC"
 export MODPATH="-module "
 export I4R4="-integer-size 32 -real-size 32"
 export I4R8="-integer-size 32 -real-size 64"
 export I8R8="-integer-size 64 -real-size 64"

 export CPPDEFS=""
 export CFLAGSDEFS=""
 export FFLAGSDEFS=""

 export USECC=""
 export USEFC="YES"
 export DEPS=""
