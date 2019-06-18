SHELL=/bin/sh
# sources, headers, objects, library and etc.
SRCS_m  = $(filter-out gdswzd_mod.f90, $(wildcard gdswzd*_mod.f90)) \
          gdswzd_mod.f90
SRCS_f  = $(SRCS_m) $(filter-out $(SRCS_m), $(wildcard *.f90))
SRCS    = $(SRCS_f) gdswzd_c.F90
OBJS    = $(patsubst %.f90, %.o, $(SRCS_f)) gdswzd_c.o

all: build

build: $(LIB)

$(LIB): $(OBJS)
	$(AR) $(ARFLAGS) $@ $^

debug: dbg_build

dbg_build: override CFLAGS += $(DEBUG)
dbg_build: override FFLAGS += $(DEBUG)
dbg_build: $(LIB)

%.o: %.c
	$(CC) -c $(CFLAGS) $(CFLAGSDEFS) $<

%.o: %.f
	$(OMPFC) -c $(FFLAGS) $(FFLAGSDEFS) $<

%.o: %.f90
	$(OMPFC) -c $(FREEFORM) $(FFLAGS) $(FFLAGSDEFS) $<

%.f90: %.F90
	$(CPP) $(CPPFLAGS) $(CPPDEFS) $< $@

message: $(MSGSRC)
	$(CC) -c $(MSGSRC)
	$(AR) $(ARFLAGS) $(LIB) $(patsubst %.c, %.o, $(MSGSRC))
	$(RM) $(MSGSRC) $(patsubst %.c, %.o, $(MSGSRC))

install:
	[[ -z "$(LIB)" || -z "$(LIB_DIR)" ]] || mv $(LIB) $(LIB_DIR)
	[[ -z "$(MOD)" || -z "$(INC_DIR)" ]] || mv $(MOD) $(INC_DIR)
	[[ -z "$(SRC_DIR)" ]] || cp ./* $(SRC_DIR)

test:
	@echo "*** unit test unavailable currently ***"

clean:
	$(RM) $(OBJS) $(LIB)
	[[ -z "$(MOD)" ]] || $(RM) -r $(MOD)

.SUFFIXES: .F90 .f90 .mod .o .a

# Specify targets that do not generate filesystem objects
.PHONY: all build install debug dbg_build message test clean
