
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

OBJ = $(SRC:.$(EXT)=.o)

CFLAGS = -I$(MODDIR) -I$(INCDIR) -fopenmp

LFLAGS  = $(SMALIB)
LFLAGS += -lpthread -lm -lgomp


ifeq ($(GPROF),yes)
   CFLAGS += -pg -g
   LFLAGS += -pg
endif


CFLAGS += -cpp -DGFORTRANMOD -ffree-form -ffree-line-length-none -march=native -fimplicit-none -fall-intrinsics -fmax-errors=1 -finit-real=nan -ffpe-summary=none
ifeq ($(DEBUG),yes)
   CFLAGS += -Og -g -Wall -Wextra -fbacktrace -pedantic -fbounds-check -Wuninitialized -fimplicit-none
else
   CFLAGS += -O3
endif


all: $(EXE)

%.o: %.$(EXT)
	$(FORT) $(CFLAGS) -c $< -o $(OBJDIR)/$@
	-mv -f *.mod $(MODDIR)

$(EXE): $(OBJ)
	$(FORT) $(addprefix $(OBJDIR)/, $(OBJ)) $(LFLAGS) -o $(EXE)

mod_data_arch.o :
mod_num_par.o : mod_data_arch.o

umfpack.o :
superlu.o :
dmumps_struc.o :
hsl_common90.o :
hsl_ddeps90.o :
hsl_ma48d.o : hsl_common90.o hsl_ddeps90.o

mod_solver.o : umfpack.o hsl_ma48d.o mod_data_arch.o superlu.o dmumps_struc.o mod_num_par.o

prg.o : mod_solver.o mod_num_par.o mod_data_arch.o

.PHONY: clean

clean:
	rm -f $(OBJDIR)/*.o
	rm -f $(MODDIR)/*.mod

