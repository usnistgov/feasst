%module shape

%{
#include "shape.h"
#ifdef FFTW_
  #include "fftw3.h"
#endif  // FFTW_
%}

%pythonnondynamic;

%include shape.h
