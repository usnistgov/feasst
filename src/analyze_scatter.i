%module analyze_scatter

%{
#include "analyze.h"
#include "analyze_scatter.h"
using namespace feasst;
%}

%pythonnondynamic;

%include analyze.h
%include analyze_scatter.h


