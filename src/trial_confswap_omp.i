%module trial_confswap_omp

%{
#include "trial.h"
#include "trial_confswap_omp.h"
using namespace feasst;
%}

%pythonnondynamic;

%include trial.h
%include trial_confswap_omp.h
