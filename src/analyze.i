%module analyze

%{
#include "analyze.h"
#include "analyze_cluster.h"
#include "analyze_orient.h"
#include "analyze_scatter.h"
#include "analyze_scatter_aniso.h"
using namespace feasst;
%}

%pythonnondynamic;

%include analyze.h
%include analyze_cluster.h
%include analyze_orient.h
%include analyze_scatter.h
%include analyze_scatter_aniso.h

