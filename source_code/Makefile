IDIR = include
CXX = g++ #-lstdc++ 
mpigXX = mpicxx
# MPI_COMPILE_FLAGS = $(shell mpicc --showme:compile)
MPI_LINK_FLAGS = $(shell mpicc --showme:link)
MPI_COMPILE_FLAGS = -fopenmp 
CXXFLAGS = -I$(IDIR) -std=c++17
ODIR = obj
_DEPS = my_functions.h particle.h beam.h threevector.h threematrix.h screen.h magnet.h
#_DEPS = threevector.h threematrix.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))
_OBJ = my_functions.o particle.o beam.o mag_spec_tracker.o screen.o magnet.o #quad_test.o #replace with mag_spec_tracker.o to test quadrupole
#_OBJ = quad_test.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

$(ODIR)/%.o: %.cpp $(DEPS)
	$(mpigXX) -c -o $@ $< $(CXXFLAGS)

run: $(OBJ)
	$(mpigXX) $(MPI_COMPILE_FLAGS) -o $@ $^ $(CXXFLAGS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(IDIR)/*~ 
