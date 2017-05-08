%module pair_wall

%{
#include "pair.h"
#include "pair_wall.h"
using namespace feasst;
%}

%pythonnondynamic;

%include pair.h
%include pair_wall.h


