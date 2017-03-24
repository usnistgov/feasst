%module trial_cluster

%{
#include "trial.h"
#include "trial_cluster.h"
%}

%pythonnondynamic;

%include trial.h
%include trial_cluster.h
