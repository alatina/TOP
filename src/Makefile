PYTHON_BIN = 
PYTHON_CONFIG = 
NUMPY_INCLUDES =  -I`$(PYTHON_BIN) -c "import numpy ; print(numpy.get_include())"`
PYTHON_CFLAGS = $(filter-out -Wstrict-prototypes,  $(shell $(PYTHON_CONFIG) --cflags)) $(NUMPY_INCLUDES)
PYTHON_LDFLAGS = $(shell $(PYTHON_CONFIG) --ldflags)
MKOCTFILE = mkoctfile

TMPFILE := $(shell mktemp)

OBJS = TOP.o

LDFLAGS = -L/opt/X11/lib -ltcl -ltk
CXXFLAGS = -I/opt/X11/include

all: TOP

install: all
	@echo "Done!"

TOP: TOP.cc 
	@$(CXX) TOP.cc -o TOP $(CXXFLAGS) $(LDFLAGS)

src*/%.o : src*/%.cc include*/%.hh
	@echo "Compiling $< ..."
	@$(CXX) $(CXXFLAGS) $(GSL_CFLAGS) $< -c -o $@

src*/%.o : src*/%.cc
	@echo "Compiling $< ..."
	@$(CXX) $(CXXFLAGS) $(GSL_CFLAGS) $< -c -o $@

src-physics/space_charge.o: include-physics/space_charge*.hh
src-physics/space_charge_pic.o: include-physics/space_charge*.hh include-math/greens_functions/*.hh

# CLEAN

clean:
	@rm -f TOP $(OBJS)\
 RF_Track_python_wrap.o\
 RF_Track_octave_wrap.o\
 RF_Track.oct\
 RF_Track.py\
 _RF_Track.so

cleanall: clean
	@rm -f RF_Track_octave_wrap.cc RF_Track_python_wrap.cc
