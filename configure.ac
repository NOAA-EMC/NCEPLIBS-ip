# Put version control id tag in output configure script
AC_REVISION($Id: configure.ac 74394 2016-04-19 19:15:36Z george.gayno@noaa.gov $)


# Initialise autoconf
AC_INIT([nceplibs ip], [v3.0.0], [ncep.list.emc.nceplibs.support@noaa.gov], [ip])


# Output info for user
AC_MSG_NOTICE([AC_PACKAGE_NAME AC_PACKAGE_VERSION])


# Check for existence of unique file before proceeding
AC_CONFIG_SRCDIR([sorc/ipolates.f90])


# Define the configuration files
AC_CONFIG_FILES([Makefile sorc/Makefile])


# Check for programs
AC_PROG_FC(ifort gfortran xlf2003 pgf95)
AC_PROG_INSTALL
AC_PROG_RANLIB
AC_CHECK_TOOL([AR], ar)
# ...Set defaults for ARFLAGS since there is no autoconf macro
test -n "$ARFLAGS" || ARFLAGS="crvs"
AC_SUBST(ARFLAGS)


# Set the programming language
AC_LANG(Fortran)
AC_FC_SRCEXT(f90)


# Establish the default type promotion for naming
AC_ARG_ENABLE([promote],
  [AS_HELP_STRING([--enable-promote],[Promote defaults to 4,8,d. @<:@default: 4@:>@])],
  [promote=${enableval}],[promote=4])
case ${promote} in
  8|d) precision=${promote};;
  *) precision=4;;
esac
AC_MSG_NOTICE([Building precision ${precision} version.])
AC_SUBST(PRECISION,[${precision}])

AC_SUBST(FPPFLAGS)

# Set platform/compiler specific flags
case $FC in
  ifort)
    case ${precision} in
      4) FPPFLAGS="-fpp -DLSIZE=4 -save-temps $FPPFLAGS";;
      8) FCFLAGS="-r8 -i8 $FCFLAGS"
         FPPFLAGS="-fpp -DLSIZE=8 -save-temps $FPPFLAGS";;
      d) FCFLAGS="-r8 $FCFLAGS"
         FPPFLAGS="-fpp -DLSIZE=d -save-temps $FPPFLAGS";;
    esac
    ;;
  pgf95)
    case ${precision} in
      4) FPPFLAGS="-DLSIZE=4 $FPPFLAGS";;
      8) FCFLAGS="-r8 -i8 $FCFLAGS"
         FPPFLAGS="-DLSIZE=8 $FPPFLAGS";;
      d) FCFLAGS="-r8 $FCFLAGS"
         FPPFLAGS="-DLSIZE=d $FPPFLAGS";;
    esac
    ;;
  gfortran)
    case ${precision} in
      4) FPPFLAGS="-cpp -DLSIZE=4 $FPPFLAGS";;
      8) FCFLAGS="-fdefault-real-8 -fdefault-integer-8 $FCFLAGS"
         FPPFLAGS="-cpp -DLSIZE=8 $FPPFLAGS";;
      d) FCFLAGS="-fdefault-real-8 $FCFLAGS"
         FPPFLAGS="-cpp -DLSIZE=d $FPPFLAGS";;
    esac
    ;;
  xlf*)
    ARFLAGS="-X64 $ARFLAGS"
    case ${precision} in
      4) FPPFLAGS="-WF,-DLSIZE=4 $FPPFLAGS";;
      8) FCFLAGS="-qrealsize=8 -qintsize=8 $FCFLAGS"
         FPPFLAGS="-WF,-DLSIZE=8 $FPPFLAGS";;
      d) FCFLAGS="-qrealsize=8 $FCFLAGS"
         FPPFLAGS="-WF,-DLSIZE=d $FPPFLAGS";;
    esac
    ;;
  *) AC_MSG_WARN([Precision promotion flags for $FC unknown and not set]);;
esac

AC_ARG_ENABLE([wcoss_cray_dir],
  [AS_HELP_STRING([--enable-wcoss_cray_dir],[Specify wcoss-cray install sub-directory. @<:@default: $FC@:>@])],
  [wcoss_cray_dir=${enableval}],[wcoss_cray_dir=$FC])
AC_SUBST(WCOSS_CRAY_DIR,[${wcoss_cray_dir}])

AC_OUTPUT
