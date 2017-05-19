%module base

%{
#include "base.h"
#include "base_random.h"
#include "base_all.h"
%}

%pythonnondynamic;

%include base.h
%include base_random.h
%include base_all.h

