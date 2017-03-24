%module trial_swap

%{
#include "trial.h"
#include "trial_swap.h"
%}

%pythonnondynamic;

%include trial.h
%include trial_swap.h
