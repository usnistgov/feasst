%module trial_grow

%{
#include "trial.h"
#include "trial_grow.h"
%}

%pythonnondynamic;

%include trial.h
%include trial_grow.h
