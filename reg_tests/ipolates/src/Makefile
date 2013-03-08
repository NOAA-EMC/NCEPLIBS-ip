# ipolates/src/Makefile.  Generated from Makefile.in by configure.

# src/ level makefile template

# Package-specific substitution variables
package = ip
version = 2.0.1
suffix  = _d

# Prefix-specific substitution variables
prefix      = /global/save/George.Gayno/ip_autoconf_regt/reg_tests

# Mixed package and prefix substitution variables
installdir  = $(prefix)/ipolates/exec

# Compiler substitution variables
FC = ifort
FCFLAGS = -r8 -check all -traceback -fpe0 -ftrapuv -assume byterecl -g -FR -openmp 
LDFLAGS = -L/global/save/George.Gayno/ip_autoconf_regt/reg_tests/lib -L/nwprod/lib 
LIBS = -lip_test_d -lsp_v2.0.1_d -lbacio_v2.0.1_4 -lw3nco_v2.0.4_d  
INSTALL         = /usr/bin/install -c
INSTALL_DATA    = ${INSTALL} -m 644
INSTALL_PROGRAM = ${INSTALL}
INSTALL_SCRIPT  = ${INSTALL}

# The file definitions. This include must occur before targets.
include make.filelist

# The targets
all: exec

exec: $(OBJ_FILES)
	$(FC) $(FCFLAGS) -o ipolates.exe $(OBJ_FILES) $(LDFLAGS) $(LIBS)

clean:
	-rm -f *.o *.exe

distclean: clean
	-rm Makefile

install:
	$(INSTALL) -d $(DESTDIR)$(installdir)
	$(INSTALL_PROGRAM) ipolates.exe $(DESTDIR)$(installdir)/ipolates$(suffix).exe

Makefile: Makefile.in ../../config.status
	cd ../.. && ./config.status ipolates/src/$@

../../config.status: ../../configure
	cd ../.. && ./config.status --recheck

# Specify targets that do not generate filesystem objects
.PHONY: all clean distclean install

# File dependency and suffix rule include files
include make.dependencies
include make.rules
