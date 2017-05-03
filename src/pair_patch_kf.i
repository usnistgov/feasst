%module pair_patch_kf

%{
#include "pair.h"
#include "pair_patch_kf.h"
using namespace feasst;
%}

%pythonnondynamic;

%include pair.h
%include pair_patch_kf.h
