%module analyze_cluster

%{
#include "analyze.h"
#include "analyze_cluster.h"
using namespace feasst;
%}

%pythonnondynamic;

%include analyze.h
%include analyze_cluster.h


