%module pair_patch_kf

%{
#include "pair.h"
#include "pair_patch_kf.h"
%}

%pythonnondynamic;

%include pair.h
%include pair_patch_kf.h
