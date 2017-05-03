%module trial_grow

%{
#include "trial.h"
#include "trial_grow.h"
using namespace feasst;
%}

%pythonnondynamic;

%include trial.h
%include trial_grow.h
