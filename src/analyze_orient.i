%module analyze_orient

%{
#include "analyze.h"
#include "analyze_orient.h"
using namespace feasst;
%}

%pythonnondynamic;

%include analyze.h
%include analyze_orient.h


