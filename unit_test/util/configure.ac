## Put version control id tag in output configure script
AC_REVISION($Id: configure.ac 25875 2013-02-19 15:45:31Z george.gayno@noaa.gov $)

## Initialise autoconf
AC_INIT([ip_unit_test], [2.0.4], [ncep.list.emc.nceplibs.support@noaa.gov])

# Output info for user
AC_MSG_NOTICE([ip AC_PACKAGE_VERSION])

# Check for existence of unique file before proceeding
# Can only check for one file.  Not a file in each
# subdirectory.
AC_CONFIG_SRCDIR([diffgb.f])

# Define the configuration files
AC_CONFIG_FILES([Makefile])

# Check for programs
AC_PROG_FC(ifort)
AC_PROG_INSTALL

# Set the programming language
AC_LANG(Fortran)
AC_FC_SRCEXT(f)

AC_OUTPUT
