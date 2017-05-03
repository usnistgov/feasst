%module trial_md

%{
#include "trial.h"
#include "trial_md.h"
using namespace feasst;
%}

%pythonnondynamic;

%include trial.h
%include trial_md.h
