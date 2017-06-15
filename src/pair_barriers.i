%module pair_barriers

%{
#include "pair.h"
#include "pair_barriers.h"
%}

%pythonnondynamic;

%include pair.h
%include pair_barriers.h
