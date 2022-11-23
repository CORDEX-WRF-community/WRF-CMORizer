FC = mpifort
FCFLAGS = -O3 -assume realloc_lhs
FCFLAGS += -fp-model precise -prec-div -prec-sqrt
FCFLAGS += -I$(EBROOTNETCDFMINFORTRAN)/include
LDFLAGS = -L$(EBROOTNETCDFMINFORTRAN)/lib -lnetcdff -lnetcdf

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
