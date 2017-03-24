%module trial_add

%{
#include "trial.h"
#include "trial_add.h"
%}

%pythonnondynamic;

%include trial.h
%include trial_add.h
