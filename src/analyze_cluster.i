%module analyze_cluster

%{
#include "analyze.h"
#include "analyze_cluster.h"
%}

%pythonnondynamic;

%include analyze.h
%include analyze_cluster.h


