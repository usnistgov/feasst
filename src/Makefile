### User settings for optional compiler flags and location of libraries  ###
#CXX = icpc
CXX = g++
CXXFLAGS = -O3 -std=c++0x -g
#CXXFLAGS += -fprofile-arcs -ftest-coverage
#CXXFLAGS = -std=c++0x -g -fprofile-arcs -ftest-coverage
#CXX = ccache g++

# # openMP
# CXXFLAGS += -DOMP_H_ -fopenmp

# # XDR library
# CXXFLAGS += -DXDRFILE_H_ -DCPLUSPLUS -I$(HOME)/include/xdrfile -L$(HOME)/lib
# LIBS += -lxdrfile

# # FFTW library
# CXXFLAGS += -DFFTW_ -I$(HOME)/software/fftw-3.3.4/build/include -L$(HOME)/software/fftw-3.3.4/build/lib 
# LIBS += -lfftw3

# # json library (requires gcc 4.9+, std=c++11 or higher)
# CXXFLAGS += -DJSON_

# # hdf5
# CXXFLAGS += -DHDF5_ -I/usr/local/hdf5/include -L/usr/local/hdf5/lib
# LIBS +=  -lhdf5 -lhdf5_cpp

# # gnu scientific library (gsl)
# CXXFLAGS += -DGSL_ -I$(HOME)/software/gsl-2.3/include -L$(HOME)/software/gsl-2.3/lib
# LIBS += -lgsl -lgslcblas -lm

# gtest
GTEST_DIR = $(HOME)/software/gtest-1.7.0
CPPFLAGS += -I$(GTEST_DIR)/include

######################################################
### User should not need to change below this line ###
######################################################

# obtain some paths
mkfile_path := $(abspath $(lastword $(MAKEFILE_LIST)))
current_dir := $(abspath $(patsubst %/,%,$(dir $(mkfile_path))))

# obtain git version
GIT_VERSION := $(shell git describe --abbrev=10 --dirty --always --tags)

CXXFLAGS += -Wall -Wextra -DFEASST_SRC_=$(current_dir) -DVERSION="\"$(GIT_VERSION)\""

base_classes = custom_exception.o shape.o base.o functions.o table.o accumulator.o histogram.o random.o space.o pair.o analyze.o criteria.o trial.o mc.o ui_abbreviated.o
SRC = $(wildcard *.cc)
OBJ = $(SRC:.cc=.o)
containing = $(foreach v,$2,$(if $(findstring $1,$v),$v))
notcontaining = $(foreach v,$2,$(if $(findstring $1,$v),,$v))
OBJnoT = $(call notcontaining,unittest,$(OBJ))
pairs = $(call containing,pair_,$(OBJnoT))
trial = $(call containing,trial_,$(OBJnoT))
crits = $(call containing,criteria_,$(OBJnoT))
anlzs = $(call containing,analyze_,$(OBJnoT))
rans = $(call containing,random_,$(OBJnoT))
bases = $(call containing,base_,$(OBJnoT))
mcs = $(call containing,mc_,$(OBJnoT))
accs = $(call containing,accumulator_,$(OBJnoT))
hists = $(call containing,histogram_,$(OBJnoT))
test_objects = $(call containing,_unittest,$(OBJ))
objects = $(base_classes) $(pairs) $(trial) $(crits) $(anlzs) $(mcs) $(rans) $(accs) $(hists) $(bases)
swigfiles = ${objects:.o=.i}
swigpyfiles = ${swigfiles:.i=.py}
GTEST_SRCS_ = $(GTEST_DIR)/src/*.cc $(GTEST_DIR)/src/*.h \
	      $(GTEST_DIR)/include/gtest/*.h \
	      $(GTEST_DIR)/include/gtest/internal/*.h

all: unittest cctest swig ui_text

c: unittest cctest

cnotest: $(objects)

main : main.o $(objects) 
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS)

ui_text : ui_text.o $(objects) 
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS)

unittest : $(objects) $(test_objects) gtest_main.a
	$(CXX) $(CXXFLAGS) -pthread -o $@ $^ $(LIBS)

gtest_main.a : gtest-all.o gtest_main.o
	$(AR) $(ARFLAGS) $@ $^	

gtest-all.o : $(GTEST_SRCS_)
	$(CXX) $(CPPFLAGS) -I$(GTEST_DIR) $(CXXFLAGS) -c $(GTEST_DIR)/src/gtest-all.cc

gtest_main.o : $(GTEST_SRCS_)
	$(CXX) $(CPPFLAGS) -I$(GTEST_DIR) $(CXXFLAGS) -c $(GTEST_DIR)/src/gtest_main.cc

#.PHONY : swignotest
swignotest : $(swigfiles)
	python setup.py build_ext --inplace "$(CXXFLAGS) $(LIBS)"

.PHONY : doxygen
doxygen :
	doxygen Doxyfile
	../tools/makelatex.sh

.PHONY : cctest
cctest : unittest
	mkdir -p tmp
#	valgrind --leak-check=full ./unittest
#	valgrind ./unittest
	./unittest --gtest_shuffle
	#./unittest 

.PHONY : swig
swig : swignotest
	./swigtest/swigtest.py

.PHONY : clean
clean: cleanSwig cleanc
	$(RM) *.o
	$(RM) -r tmp
	$(RM) .depend
	$(RM) -r ../html ../latex
	
.PHONY : cleanc
cleanc:
	$(RM) *.o
	$(RM) gtest_main.a gtest-all.o gtest_main.o
	$(RM) main
	$(RM) unittest
	$(RM) ui_text
	$(RM) gmon.out unittest
	$(RM) *.gcda *.gcno *.gcov

.PHONY : cleanSwig
cleanSwig:
	$(RM) -r build
	$(RM) *_wrap.cxx
	$(RM) *_wrap.cpp
	$(RM) feasst.py
	$(RM) feasst.pyc
	$(RM) pyfeasst.pyc
	$(RM) _feasst.so
	$(RM) $(swigpyfiles)	


depend:
	$(CXX) -MM $(CXXFLAGS) *.cc > .depend

-include .depend


