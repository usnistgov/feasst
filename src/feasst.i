%module feasst

%ignore CriteriaWLTMMC::lnPIrwsatwrap;
%ignore CriteriaWLTMMC::lnPIrwnmxwrap;
%ignore MC::boyleminwrap;

%{
#include "functions.h"
#include "custom_exception.h"
#include "histogram.h"
#include "accumulator.h"
#include "accumulator_vec.h"
#include "analyze.h"
#include "analyze_cluster.h"
#include "analyze_scatter.h"
#include "base.h"
#include "base_random.h"
#include "base_all.h"
#include "space.h"
#ifdef XDRFILE_H_
  extern "C" {
    #include "xdrfile.h"
    #include "xdrfile_xtc.h"
    #include "xdrfile_trr.h"
  }
#endif // XDRFILE_H_
#include "pair.h"
#include "pair_lj.h"
#include "pair_lj_multi.h"
#include "pair_patch_kf.h"
#include "pair_lj_coul_ewald.h"
#include "pair_ideal.h"
#include "pair_hard_circle.h"
#include "pair_squarewell.h"
#include "pair_tabular_1d.h"
#include "pair_hybrid.h"
#include "pair_hs.h"
#include "trial.h"
#include "trial_add.h"
#include "trial_delete.h"
#include "trial_transform.h"
#include "trial_confswap_omp.h"
#include "trial_confswap_txt.h"
#include "trial_xswap.h"
#include "trial_md.h"
#include "trial_swap.h"
#include "criteria.h"
#include "criteria_mayer.h"
#include "criteria_metropolis.h"
#include "criteria_wltmmc.h"
#include "mc.h"
#include "mc_wltmmc.h"
#include "ui_abbreviated.h"
%}

%include "std_vector.i"
%template(IntVector) std::vector<int>;
%template(DoubleVector) std::vector<double>;
using namespace std;
// using std::vector;

%pythonnondynamic;

%include functions.h
%include custom_exception.h
%include histogram.h
%include accumulator.h
%include accumulator_vec.h
%include analyze.h
%include analyze_cluster.h
%include analyze_scatter.h
%include base.h
%include base_random.h
%include base_all.h
%include space.h
%include pair.h
%include pair_lj.h
%include pair_lj_multi.h
%include pair_patch_kf.h
%include pair_lj_coul_ewald.h
%include pair_ideal.h
%include pair_hard_circle.h
%include pair_squarewell.h
%include pair_tabular_1d.h
%include pair_hybrid.h
%include pair_hs.h
%include trial.h
%include trial_add.h
%include trial_delete.h
%include trial_transform.h
%include trial_confswap_omp.h
%include trial_confswap_txt.h
%include trial_xswap.h
%include trial_md.h
%include trial_swap.h
%include criteria.h
%include criteria_mayer.h
%include criteria_metropolis.h
%include criteria_wltmmc.h
%include mc.h
%include mc_wltmmc.h
%include ui_abbreviated.h



