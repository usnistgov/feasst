%module trial_transform

%{
#include "trial.h"
#include "trial_transform.h"
using namespace feasst;
%}

%pythonnondynamic;

%include trial.h
%include trial_transform.h
