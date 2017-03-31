%module space

%{
#include "space.h"
#ifdef XDRFILE_H_
  extern "C" {
    #include "xdrfile.h"
    #include "xdrfile_xtc.h"
    #include "xdrfile_trr.h"
  }
#endif  // XDRFILE_H_
%}

%pythonnondynamic;

%include space.h
