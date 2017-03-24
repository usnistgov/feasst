%module pair_lj

%{
#include "pair.h"
#include "pair_lj.h"
#include "pair_lj_multi.h"
%}

%pythonnondynamic;

%include pair.h
%include pair_lj.h
%include pair_lj_multi.h


