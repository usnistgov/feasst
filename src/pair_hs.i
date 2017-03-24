%module pair_hs

%{
#include "pair.h"
#include "pair_hs.h"
%}

%pythonnondynamic;

%include pair.h
%include pair_hs.h


