%module barrier

%{
#include "barrier.h"
#include "barrier_planar.h"
%}

%pythonnondynamic;

%include barrier.h
%include barrier_planar.h

