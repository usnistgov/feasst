%module trial_pressure

%{
#include "trial.h"
#include "trial_pressure.h"
using namespace feasst;
%}

%pythonnondynamic;

%include trial.h
%include trial_pressure.h
