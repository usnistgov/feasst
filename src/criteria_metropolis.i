%module criteria_metropolis

%{
#include "criteria.h"
#include "criteria_metropolis.h"
using namespace feasst;
%}

%pythonnondynamic;

%include criteria.h
%include criteria_metropolis.h
