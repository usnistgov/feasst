%module analyze_scatter_aniso

%{
#include "analyze.h"
#include "analyze_scatter_aniso.h"
#ifdef FFTW_
  #include "fftw3.h"
#endif  // FFTW_
using namespace feasst;
%}

%pythonnondynamic;

%include analyze.h
%include analyze_scatter_aniso.h


