%module accumulator

%{
#include "accumulator.h"
#include "accumulator_vec.h"
%}

%pythonnondynamic;

%include accumulator.h
%include accumulator_vec.h

