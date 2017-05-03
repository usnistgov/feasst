%module trial_add

%{
#include "trial.h"
#include "trial_add.h"
using namespace feasst;
%}

%pythonnondynamic;

%include trial.h
%include trial_add.h
