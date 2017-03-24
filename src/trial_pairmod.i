%module trial_pairmod

%{
#include "trial.h"
#include "trial_pairmod.h"
%}

%pythonnondynamic;

%include trial.h
%include trial_pairmod.h
