## Put version control id tag in output configure script
AC_REVISION($Id: configure.ac 25875 2013-02-19 15:45:31Z george.gayno@noaa.gov $)

## Initialise autoconf
AC_INIT([ip regression test suite], [2.0.4], [ncep.list.emc.nceplibs.support@noaa.gov])

# Output info for user
AC_MSG_NOTICE([ip AC_PACKAGE_VERSION])

# Check for existence of unique file before proceeding
# Can only check for one file.  Not a file in each
# subdirectory.
AC_CONFIG_SRCDIR([gausslat/sorc/gausslat_driver.f90])

# Define the configuration files
AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([ipxetas/sorc/Makefile])
AC_CONFIG_FILES([ipxwafs/sorc/Makefile])
AC_CONFIG_FILES([copygb2/sorc/Makefile])
AC_CONFIG_FILES([gausslat/sorc/Makefile])
AC_CONFIG_FILES([gdswzd/sorc/Makefile])
AC_CONFIG_FILES([ipolates/sorc/Makefile])
AC_CONFIG_FILES([ipolatev/sorc/Makefile])

# Check for programs
AC_PROG_FC(ifort)
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
  ifort)
    case ${precision} in
      4) ;;
      8) FCFLAGS="-r8 -i8 $FCFLAGS";;
      d) FCFLAGS="-r8 $FCFLAGS";;
    esac
    ;;
  *) AC_MSG_WARN([Precision promotion flags for $FC unknown and not set]);;
esac

AC_OUTPUT
