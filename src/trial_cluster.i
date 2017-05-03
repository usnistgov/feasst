%module trial_cluster

%{
#include "trial.h"
#include "trial_cluster.h"
using namespace feasst;
%}

%pythonnondynamic;

%include trial.h
%include trial_cluster.h
