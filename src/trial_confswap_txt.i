%module trial_confswap_txt

%{
#include "trial.h"
#include "trial_confswap_txt.h"
using namespace feasst;
%}

%pythonnondynamic;

%include trial.h
%include trial_confswap_txt.h
