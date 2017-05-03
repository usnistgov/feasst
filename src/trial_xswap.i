%module trial_xswap

%{
#include "trial.h"
#include "trial_xswap.h"
using namespace feasst;
%}

%pythonnondynamic;

%include trial.h
%include trial_xswap.h
