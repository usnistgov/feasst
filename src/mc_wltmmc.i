%module mc_wltmmc

%ignore feasst::MC::boyleminwrap;

%{
#include "mc.h"
#include "mc_wltmmc.h"
using namespace feasst;
%}

%pythonnondynamic;

%include mc.h
%include mc_wltmmc.h
