%module base_all

%{
#include "base.h"
#include "base_random.h"
#include "base_math.h"
#include "base_all.h"
%}

%pythonnondynamic;

%include base.h
%include base_random.h
%include base_math.h
%include base_all.h


