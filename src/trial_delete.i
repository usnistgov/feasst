%module trial_delete

%{
#include "trial.h"
#include "trial_delete.h"
using namespace feasst;
%}

%pythonnondynamic;

%include trial.h
%include trial_delete.h
