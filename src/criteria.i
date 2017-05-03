%module criteria

%ignore feasst::CriteriaWLTMMC::lnPIrwsatwrap;
%ignore feasst::CriteriaWLTMMC::lnPIrwnmxwrap;

%{
#include "criteria.h"
#include "criteria_metropolis.h"
#include "criteria_wltmmc.h"
using namespace feasst;
%}

%pythonnondynamic;

%include criteria.h
%include criteria_metropolis.h
%include criteria_wltmmc.h
