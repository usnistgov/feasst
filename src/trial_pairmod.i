%module trial_pairmod

%{
#include "trial.h"
#include "trial_pairmod.h"
using namespace feasst;
%}

%pythonnondynamic;

%include trial.h
%include trial_pairmod.h
