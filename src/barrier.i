%module barrier

%{
#include "barrier.h"
#include "pair.h"
%}

%pythonnondynamic;

%include barrier.h
%include pair.h
%include base.h
