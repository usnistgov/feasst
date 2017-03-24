### User settings for optional compiler flags and location of libraries  ###
#CXX = icpc
HOST = $(shell hostname)
CXX = g++ -O3 -std=c++11 -g
ifeq ($(HOST), raritan-2)
CXX = g++ -O3 -std=c++0x
#CXX = g++ -O3 -pg -g -std=c++0x
endif
ifeq ($(HOST), pn101924.nist.gov)
#CXX = g++ -O3 -std=c++11 -g
#CXX = g++ -03 -fprofile-arcs -ftest-coverage -g -std=c++11
#CXX = g++ -fprofile-arcs -ftest-coverage -g -std=c++11
#CXX = g++ -O3 -std=c++0x
CXX = ccache g++ -O3 -pg -g -std=c++11
endif
ifeq ($(HOST), gibbs.nist.gov)
CXX = g++ -O3 -std=c++0x
endif
ifeq ($(HOST), raritan.nist.gov)
CXX = g++ -O3 -std=c++0x
endif
GTEST_DIR = $(HOME)/software/gtest-1.7.0

# obtain some paths
mkfile_path := $(abspath $(lastword $(MAKEFILE_LIST)))
current_dir := $(abspath $(patsubst %/,%,$(dir $(mkfile_path))))
#current_dir := $(notdir $(patsubst %/,%,$(dir $(mkfile_path))))

# obtain git version
GIT_VERSION := $(shell git describe --abbrev=10 --dirty --always --tags)

CXXFLAGS += -Wall -Wextra -DFEASST_SRC_=$(current_dir) -DVERSION="\"$(GIT_VERSION)\""

# openMP
CXXFLAGS += -DOMP_H_ -fopenmp

# XDR library
CXXFLAGS += -DXDRFILE_H_ -DCPLUSPLUS -I$(HOME)/include/xdrfile -L$(HOME)/lib
LIBS += -lxdrfile

# FFTW library
#CXXFLAGS += -DFFTW_ -I$(HOME)/software/fftw-3.3.4/build/include -L$(HOME)/software/fftw-3.3.4/build/lib 
#LIBS += -lfftw3

# json library (requires gcc 4.9+, std=c++11 or higher)
#CXXFLAGS += -DJSON_

# hdf5
#CXXFLAGS += -DHDF5_ -I/usr/local/hdf5/include -L/usr/local/hdf5/lib
#LIBS +=  -lhdf5 -lhdf5_cpp

# gtest
CPPFLAGS += -I$(GTEST_DIR)/include
### NOTE: SWIG requires setup.py to compile correctly ###

### User should not change below this line ###
###  Except: you can change the order of the tests by changing the order of the base_classes below ###
base_classes = custom_exception.o shape.o base.o functions.o table.o accumulator.o histogram.o random.o space.o pair.o analyze.o criteria.o trial.o mc.o ui_abbreviated.o
#base_classes = mc.o trial.o criteria.o analyze.o pair.o space.o random.o histogram.o accumulator.o table.o functions.o base.o shape.o
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
#test_objects = ${base_classes:.o=_unittest.o}
objects = $(base_classes) $(pairs) $(trial) $(crits) $(anlzs) $(mcs) $(rans) $(accs) $(hists) $(bases)
swigfiles = ${objects:.o=.i}
swigpyfiles = ${swigfiles:.i=.py}
GTEST_SRCS_ = $(GTEST_DIR)/src/*.cc $(GTEST_DIR)/src/*.h \
	      $(GTEST_DIR)/include/gtest/*.h \
	      $(GTEST_DIR)/include/gtest/internal/*.h

all: unittest doxygen swig cctest swigtest

c: unittest cctest

cnotest: $(objects)

main : main.o $(objects) 
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS)

unittest : $(objects) $(test_objects) gtest_main.a
	$(CXX) $(CXXFLAGS) -pthread -o $@ $^ $(LIBS)

gtest_main.a : gtest-all.o gtest_main.o
	$(AR) $(ARFLAGS) $@ $^	

gtest-all.o : $(GTEST_SRCS_)
	$(CXX) $(CPPFLAGS) -I$(GTEST_DIR) $(CXXFLAGS) -c $(GTEST_DIR)/src/gtest-all.cc

gtest_main.o : $(GTEST_SRCS_)
	$(CXX) $(CPPFLAGS) -I$(GTEST_DIR) $(CXXFLAGS) -c $(GTEST_DIR)/src/gtest_main.cc

#.PHONY : swig
swig : $(swigfiles)
	#swig -3ython -c++ feasst.i
	python setup.py build_ext --inplace $(current_dir) "\"$(GIT_VERSION)\""
#%.i : %.o
#	swig -python -c++ $@

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

.PHONY : swigtest
swigtest : swig
	./swigtest.py

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
	$(RM) gmon.out unittest
	$(RM) *.gcda *.gcno

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


