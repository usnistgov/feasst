%module trial_swap

%{
#include "trial.h"
#include "trial_swap.h"
using namespace feasst;
%}

%pythonnondynamic;

%include trial.h
%include trial_swap.h
