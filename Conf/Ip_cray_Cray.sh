# *** for WCOSS Cray (cray) ***

 export CC=cc
 export FC=ftn
 export CPP=cpp
 export OMPCC="$CC -homp"
 export OMPFC="$FC -homp"
 export MPICC=$CC
 export MPIFC=$FC

 export DEBUG=""
 export CFLAGS=""
 export FFLAGS="-O2 -homp -hoverindex -hPIC -I."
 export FPPCPP=""
 export FREEFORM=""
 export CPPFLAGS=""
 export MPICFLAGS=$CFLAGS
 export MPIFFLAGS=$FFLAGS
 export MODPATH="-J "
 export I4R4="-s integer32 -s real32"
 export I4R8="-s integer32 -s real64"
 export I8R8="-s integer64 -s real64"

 export CPPDEFS=""
 export CFLAGSDEFS="-DUNDERSCORE -DLINUX"
 export FFLAGSDEFS=""

 export USECC=""
 export USEFC="YES"
 export DEPS=""

 export VERSION_FLAG="-V"
