%module trial_avb

%{
#include "trial.h"
#include "trial_avb.h"
using namespace feasst;
%}

%pythonnondynamic;

%include trial.h
%include trial_avb.h
