%module criteria

%ignore CriteriaWLTMMC::lnPIrwsatwrap;
%ignore CriteriaWLTMMC::lnPIrwnmxwrap;

%{
#include "criteria.h"
#include "criteria_metropolis.h"
#include "criteria_wltmmc.h"
#include "criteria_mayer.h"
%}

%pythonnondynamic;

%include criteria.h
%include criteria_metropolis.h
%include criteria_wltmmc.h
%include criteria_mayer.h
