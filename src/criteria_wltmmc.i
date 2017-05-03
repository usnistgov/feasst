%module criteria_wltmmc

%ignore feasst::CriteriaWLTMMC::lnPIrwsatwrap;
%ignore feasst::CriteriaWLTMMC::lnPIrwnmxwrap;

%{
#include "criteria.h"
#include "criteria_wltmmc.h"
using namespace feasst;
%}

%pythonnondynamic;

%include criteria.h
%include criteria_wltmmc.h
