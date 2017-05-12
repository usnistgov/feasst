%module base

%{
#include "base.h"
#include "base_random.h"
#include "base_all.h"
using namespace feasst;
%}

%pythonnondynamic;

%include base.h
%include base_random.h
%include base_all.h

