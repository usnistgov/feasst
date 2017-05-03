%module pair_hybrid

%{
#include "pair.h"
#include "pair_hybrid.h"
using namespace feasst;
%}

%pythonnondynamic;

%include pair.h
%include pair_hybrid.h


