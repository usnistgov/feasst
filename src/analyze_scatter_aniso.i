%module analyze_scatter_aniso

%{
#include "analyze.h"
#include "analyze_scatter_aniso.h"
#include "fftw3.h"
%}

%pythonnondynamic;

%include analyze.h
%include analyze_scatter_aniso.h


