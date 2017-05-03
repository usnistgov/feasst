%module pair_lj_multi

%{
#include "pair.h"
#include "pair_lj.h"
#include "pair_lj_multi.h"
using namespace feasst;
%}

%pythonnondynamic;

%include pair.h
%include pair_lj.h
%include pair_lj_multi.h


