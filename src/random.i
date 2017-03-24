%module random

%{
#include "random.h"
#include "random_nr3.h"
%}

%pythonnondynamic;

%include random.h
%include random_nr3.h

