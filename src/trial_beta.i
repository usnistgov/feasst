%module trial_beta

%{
#include "trial.h"
#include "trial_beta.h"
using namespace feasst;
%}

%pythonnondynamic;

%include trial.h
%include trial_beta.h
