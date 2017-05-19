%module trial

%{
#include "trial.h"
#include "trial_add.h"
#include "trial_delete.h"
#include "trial_avb.h"
#include "trial_transform.h"
#include "trial_configBias.h"
#include "trial_confswap_omp.h"
#include "trial_confswap_txt.h"
#include "trial_grow.h"
#include "trial_gca.h"
#include "trial_pairmod.h"
#include "trial_beta.h"
#include "trial_cluster.h"
#include "trial_pressure.h"
#include "trial_xswap.h"
#include "trial_md.h"
#include "trial_swap.h"
%}

%pythonnondynamic;

%include trial.h
%include trial_add.h
%include trial_delete.h
%include trial_avb.h
%include trial_transform.h
%include trial_configBias.h
%include trial_confswap_omp.h
%include trial_confswap_txt.h
%include trial_grow.h
%include trial_gca.h
%include trial_pairmod.h
%include trial_beta.h
%include trial_cluster.h
%include trial_pressure.h
%include trial_xswap.h
%include trial_md.h
%include trial_swap.h
