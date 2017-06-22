/**
 * Factory method for Pair, Criteria, Trial and Analyze.
 * This implementation resides in a special file so that it may be changed
 * and diverge for different branches/repos with different files.
 */

#include "./pair_lj_multi.h"
#include "./pair_lj_coul_ewald.h"
#include "./pair_tabular_1d.h"
#include "./pair_hs.h"

#include "./criteria_metropolis.h"
#include "./criteria_wltmmc.h"
#include "./criteria_mayer.h"

#include "./trial_transform.h"
#include "./trial_add.h"
#include "./trial_delete.h"
#include "./trial_xswap.h"
#ifdef MPI_H_
  #include "./trial_confswap_txt.h"
#endif  // MPI_H_
#ifdef _OPENMP
  #include "./trial_confswap_omp.h"
#endif  // _OPENMP

#include "./analyze.h"

#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

Pair* Pair::makePair(Space* space,
  const char* fileName
  ) {
  ASSERT(fileExists(fileName),
    "restart file(" << fileName << ") doesn't exist");
  string pairtypestr = fstos("className", fileName);
  Pair* pair = NULL;
  if (pairtypestr.compare("PairLJ") == 0) {
    pair = new PairLJ(space, fileName);
  } else if (pairtypestr.compare("PairLJMulti") == 0) {
    pair = new PairLJMulti(space, fileName);
  } else if (pairtypestr.compare("PairLJCoulEwald") == 0) {
    pair = new PairLJCoulEwald(space, fileName);
  } else if (pairtypestr.compare("PairTabular1D") == 0) {
    pair = new PairTabular1D(space, fileName);
  } else if (pairtypestr.compare("PairHS") == 0) {
    pair = new PairHS(space, fileName);
//  } else if (pairtypestr.compare("PairPatchHPC") == 0) {
//    pair = new PairPatchHPC(space, fileName);
//  } else if (pairtypestr.compare("PairPatchKF") == 0) {
//    pair = new PairPatchKF(space, fileName);
//  } else if (pairtypestr.compare("PairIdeal") == 0) {
//    pair = new PairIdeal(space, fileName);
  } else {
    ASSERT(0, "unrecognized pair(" << pairtypestr << ") in factory");
  }
  return pair;
}

Criteria* Criteria::makeCriteria(const char* fileName
  ) {
  ASSERT(fileExists(fileName),
    "restart file(" << fileName << ") doesn't exist");
  string typestr = fstos("className", fileName);
  Criteria* criteria = NULL;
  if ( (typestr.compare("Criteria") == 0) ||
       (typestr.compare("CriteriaMetropolis") == 0) ) {
    criteria = new CriteriaMetropolis(fileName);
  } else if (typestr.compare("CriteriaWLTMMC") == 0) {
    criteria = new CriteriaWLTMMC(fileName);
  } else if (typestr.compare("CriteriaMayer") == 0) {
    criteria = new CriteriaMayer(fileName);
  } else {
    ASSERT(0, "unrecognized criteria(" << typestr << ") in factory");
  }
  return criteria;
}

shared_ptr<Trial> Trial::makeTrial(
  Space* space,
  Pair* pair,
  Criteria* criteria,
  const char* fileName) {
  ASSERT(fileExists(fileName),
         "restart file(" << fileName << ") doesn't exist");
  string trialtypestr = fstos("className", fileName);
  shared_ptr<Trial> trial;
  if (trialtypestr.compare("TrialTransform") == 0) {
    trial = make_shared<TrialTransform>(fileName, space, pair, criteria);
  } else if (trialtypestr.compare("TrialAdd") == 0) {
    trial = make_shared<TrialAdd>(fileName, space, pair, criteria);
  } else if (trialtypestr.compare("TrialDelete") == 0) {
    trial = make_shared<TrialDelete>(fileName, space, pair, criteria);
#ifdef _OPENMP
  } else if (trialtypestr.compare("TrialConfSwapOMP") == 0) {
    trial = make_shared<TrialConfSwapOMP>(fileName, space, pair, criteria);
#endif  // _OPENMP
#ifdef MPI_H_
  } else if (trialtypestr.compare("TrialConfSwapTXT") == 0) {
    trial = make_shared<TrialConfSwapTXT>(fileName, space, pair, criteria);
#endif  // MPI_H_
  } else if (trialtypestr.compare("TrialXSwap") == 0) {
    trial = make_shared<TrialXSwap>(fileName, space, pair, criteria);
  } else {
    ASSERT(0, "unrecognized trial class(" << trialtypestr);
  }
  return trial;
}

shared_ptr<Analyze> Analyze::makeAnalyze(
  Space* space,
  Pair* pair,
  const char* fileName) {
  ASSERT(fileExists(fileName),
         "restart file(" << fileName << ") doesn't exist");
  string anatypestr = fstos("className", fileName);
  shared_ptr<Analyze> ana;
  ASSERT(0, "unrecognized analyze class(" << anatypestr << ")");
  return ana;
}

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_

