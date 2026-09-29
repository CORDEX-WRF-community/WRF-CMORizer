ROOT_NETCDF_FORTRAN="full-path-to-netcdf-fortran library"
ROOT_NETCDF_C="full-path-to-netcdf-c library"

# gfortran flags (uncomment if using gfortran)
#FC = gfortran
#FCFLAGS = -O2 -cpp -DSERIAL
#FCFLAGS += -Wall
#FCFLAGS += -ffree-line-length-none
#FCFLAGS += -Wno-tabs -Wno-unused-variable -Wno-maybe-uninitialized

# intel flags (comment out if not using intel compiler)
FC = ifort 
FCFLAGS = -O2 -assume realloc_lhs -cpp -DSERIAL
FCFLAGS += -fp-model precise -prec-div -prec-sqrt

# Flags for gfortran and intel (never comment out)
FCFLAGS += -I$(ROOT_NETCDF_FORTRAN)/include -I$(ROOT_NETCDF_C)/include
LDFLAGS = -L$(ROOT_NETCDF_FORTRAN)/lib -lnetcdff -L${ROOT_NETCDF_C}/lib -lnetcdf

PROGRAMS = pCMORizer

all: $(PROGRAMS)

%: %.o
	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.f90
	$(FC) $(FCFLAGS) -c $<

.PHONY: clean veryclean

clean:
	rm -f *.o *.mod *.MOD *_genmod.f90

veryclean: clean
	rm -rf *~ $(PROGRAMS)

