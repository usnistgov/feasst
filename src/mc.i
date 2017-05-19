%module mc

%ignore MC::boyleminwrap;

%{
#include "mc.h"
#include "mc_wltmmc.h"
%}

%pythonnondynamic;

%include mc.h
%include mc_wltmmc.h
