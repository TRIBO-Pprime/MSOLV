.PHONY: clean debug prod gprof all

clean:
	rm -f $(OBJDIR)/*.o
	rm -f $(MODDIR)/*.mod
	rm prg

common_dir = ../common

MODDIR = mod
OBJDIR = obj
SRCDIR = src
LIBDIR = lib
INCDIR = inc

GENLIB  =
GENLIB += $(LIBDIR)/librefblas.a
GENLIB += $(LIBDIR)/liblapack.a

MUMLIB  =
MUMLIB += $(LIBDIR)/libmpiseq.a
MUMLIB += $(LIBDIR)/libdmumps.a
MUMLIB += $(LIBDIR)/libmumps_common.a
MUMLIB += $(LIBDIR)/libmetis.a

MUMLIB += $(LIBDIR)/libscotch.a
MUMLIB += $(LIBDIR)/libpord.a
MUMLIB += $(LIBDIR)/libesmumps.a
MUMLIB += $(LIBDIR)/libscotcherr.a

UMFLIB  =
UMFLIB += $(LIBDIR)/libumfpack.a
UMFLIB += $(LIBDIR)/libamd.a
UMFLIB += $(LIBDIR)/libsuitesparseconfig.a
UMFLIB += $(LIBDIR)/libcholmod.a
UMFLIB += $(LIBDIR)/libcolamd.a
UMFLIB += $(LIBDIR)/libccolamd.a
UMFLIB += $(LIBDIR)/libcamd.a

SUPLIB  =
SUPLIB += $(LIBDIR)/libsuperlu.a

M48LIB  =
M48LIB += $(LIBDIR)/libhsl_ma48.a

SMALIB = -Wl,--start-group $(GENLIB) $(M48LIB) $(MUMLIB) $(SUPLIB) $(UMFLIB) $(FFTLIB) -Wl,--end-group

EXT = f90

VPATH = $(MODDIR):$(OBJDIR):$(SRCDIR)

EXE = prg

SRC = $(notdir $(wildcard $(SRCDIR)/*.$(EXT)))

OBJ      = $(SRC:.$(EXT)=.o)
ALL_OBJS = $(addprefix $(OBJDIR)/, *.o) $(addprefix $(common_dir)/$(OBJDIR)/, *.o) 


CFLAGS  = -I$(MODDIR) -I$(INCDIR) -I$(common_dir)/$(MODDIR)
CFLAGS += -cpp -DGFORTRANMOD -ffree-form -ffree-line-length-none -march=native -fimplicit-none -fall-intrinsics -fmax-errors=1 -finit-real=nan -ffpe-summary=none -fopenmp

LFLAGS  = $(SMALIB)
LFLAGS += -lpthread -lm -lgomp


ifneq ('$(DEBUG)','')
	CFLAGS += -Og -g -Wall -Wextra -fbacktrace -pedantic -fbounds-check -Wuninitialized -fimplicit-none
else
	CFLAGS += -O3
endif

ifneq ('$(GPROF)','')
	CFLAGS += -pg -g
	LFLAGS += -pg
endif

all: $(EXE)

%.o: %.$(EXT)
	$(FORT) $(CFLAGS) -c $< -o $(OBJDIR)/$@
	-mv -f *.mod $(MODDIR)

$(EXE): $(OBJ)
	$(FORT) $(ALL_OBJS) $(LFLAGS) -o $(EXE)

umfpack.o :
superlu.o :
dmumps_struc.o :
hsl_common90.o :
hsl_ddeps90.o :
hsl_ma48d.o : hsl_common90.o hsl_ddeps90.o

mod_solver.o : umfpack.o hsl_ma48d.o superlu.o dmumps_struc.o

prg.o : mod_solver.o

