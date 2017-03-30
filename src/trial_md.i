%module trial_md

%{
#include "trial.h"
#include "trial_md.h"
%}

%pythonnondynamic;

%include trial.h
%include trial_md.h
