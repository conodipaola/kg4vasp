F90=gfortran
#F90=ifort
DEPS = 
FLAGS=
OBJ = preclib.o vars.o read_input.o read_dat.o sigma_real.o sum_rule.o sigma_real_dc.o kg4vasp.o

%.o: %.f90 $(DEPS)
	$(F90) -c -o $@ $< $(FLAGS)

kg4vasp.x: $(OBJ)
	$(F90) -o $@ $^ $(FLAGS)


PHONY: clean

clean:
	rm -f *.o *.mod *.x
