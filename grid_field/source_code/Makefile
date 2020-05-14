IDIR = include
CXX = g++ 
CXXFLAGS = -I$(IDIR) -std=c++17
ODIR = obj
_DEPS = my_functions.h particle.h threevector.h threematrix.h screen.h magnet.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))
_OBJ = mag_spec_tracker.o particle.o screen.o magnet.o my_functions.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

$(ODIR)/%.o: %.cpp $(DEPS)
	$(CXX) -c -o $@ $< $(CXXFLAGS)

run: $(OBJ)
	$(CXX) -o $@ $^ $(CXXFLAGS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(IDIR)/*~ 
