%module trial_beta

%{
#include "trial.h"
#include "trial_beta.h"
%}

%pythonnondynamic;

%include trial.h
%include trial_beta.h
