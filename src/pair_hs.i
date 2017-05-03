%module pair_hs

%{
#include "pair.h"
#include "pair_hs.h"
using namespace feasst;
%}

%pythonnondynamic;

%include pair.h
%include pair_hs.h


