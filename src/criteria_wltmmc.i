%module criteria_wltmmc

%ignore CriteriaWLTMMC::lnPIrwsatwrap;
%ignore CriteriaWLTMMC::lnPIrwnmxwrap;

%{
#include "criteria.h"
#include "criteria_wltmmc.h"
%}

%pythonnondynamic;

%include criteria.h
%include criteria_wltmmc.h
