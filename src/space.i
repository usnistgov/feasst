%module space

%{
#include "space.h"
extern "C" {
#include "xdrfile.h"
#include "xdrfile_xtc.h"
#include "xdrfile_trr.h"
}
%}

%pythonnondynamic;

%include space.h
