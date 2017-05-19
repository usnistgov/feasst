%module trial_pressure

%{
#include "trial.h"
#include "trial_pressure.h"
%}

%pythonnondynamic;

%include trial.h
%include trial_pressure.h
