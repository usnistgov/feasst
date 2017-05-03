%module pair_lj_coul_ewald

%{
#include "pair.h"
#include "pair_lj_coul_ewald.h"
using namespace feasst;
%}

%pythonnondynamic;

%include pair.h
%include pair_lj_coul_ewald.h
