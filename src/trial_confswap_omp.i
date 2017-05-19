%module trial_confswap_omp

%{
#include "trial.h"
#include "trial_confswap_omp.h"
%}

%pythonnondynamic;

%include trial.h
%include trial_confswap_omp.h
