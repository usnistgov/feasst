%module trial_transform

%{
#include "trial.h"
#include "trial_transform.h"
%}

%pythonnondynamic;

%include trial.h
%include trial_transform.h
