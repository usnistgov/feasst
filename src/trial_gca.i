%module trial_gca

%{
#include "trial.h"
#include "trial_gca.h"
using namespace feasst;
%}

%pythonnondynamic;

%include trial.h
%include trial_gca.h
