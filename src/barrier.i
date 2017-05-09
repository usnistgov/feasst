%module barrier

%{
#include "barrier.h"
#include "barrier_planar.h"
using namespace feasst;
%}

%pythonnondynamic;

%include barrier.h
%include barrier_planar.h

