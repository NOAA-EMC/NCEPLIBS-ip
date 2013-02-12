## Put version control id tag in output configure script
AC_REVISION($Id: configure.ac 24354 2013-01-10 20:04:24Z george.gayno@noaa.gov $)

## Initialise autoconf
AC_INIT([ip], [2.0.1], [ncep.list.emc.nceplibs.support@noaa.gov])

# Output info for user
AC_MSG_NOTICE([ip AC_PACKAGE_VERSION])

# Check for existence of unique file before proceeding
AC_CONFIG_SRCDIR([src/ipolates.f])

# Define the configuration files
AC_CONFIG_FILES([Makefile src/Makefile])

# Check for programs
AC_PROG_FC(ifort gfortran xlf2003 pgf95)
AC_PROG_INSTALL

# Set the programming language
AC_LANG(Fortran)
AC_FC_SRCEXT(f)

# Establish the default type promotion for naming
AC_ARG_ENABLE([promote],
  [AS_HELP_STRING([--enable-promote],[Promote defaults to 4,8,d. @<:@default: 4@:>@])],
  [promote=${enableval}],[promote=4])
case ${promote} in
  8|d) precision=${promote};;
  *) precision=4;;
esac
AC_MSG_NOTICE([Building precision ${precision} version.])
AC_SUBST(SUFFIX,[_${precision}])

# Set the default type promotion compiler flags
case $FC in
  ifort|gfortran|pgf95)
    case ${precision} in
      4) ;;
      8) FCFLAGS="-r8 -i8 $FCFLAGS";;
      d) FCFLAGS="-r8 $FCFLAGS";;
    esac
    ;;
  xlf2003)
    case ${precision} in
      4) ;;
      8) FCFLAGS="-qrealsize=8 -qintsize=8 $FCFLAGS";;
      d) FCFLAGS="-qrealsize=8 $FCFLAGS";;
    esac
    ;;
  *) AC_MSG_WARN([Precision promotion flags for $FC unknown and not set]);;
esac

AC_OUTPUT
