%module custom_exception

%{
#include "custom_exception.h"
using namespace feasst;
%}

%pythonnondynamic;

%include custom_exception.h

