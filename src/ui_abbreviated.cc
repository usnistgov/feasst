#ifndef UI_ABBREVIATED_H_
#define UI_ABBREVIATED_H_

#include "./ui_abbreviated.h"
#include "./trial_transform.h"
#include "./trial_add.h"
#include "./trial_delete.h"
#include "./trial_avb.h"
#include "./trial_grow.h"
#include "./trial_configBias.h"
#include "./trial_gca.h"
#include "./trial_pairmod.h"
#include "./trial_beta.h"
#include "./trial_cluster.h"
// #include "./trial_pressure.h"
#include "./trial_xswap.h"
#include "./trial_swap.h"
#include "./mc.h"
#include "./mc_wltmmc.h"

//using namespace feasst;
#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

void transformTrial(MC *mc, const char* type, double maxMoveParam = -1) {
  string typeStr(type);

  // the following if statements are flags to catch legacy code
  if (typeStr.compare("gca") == 0) {
    shared_ptr<TrialGCA> trial = make_shared<TrialGCA>();
    if (maxMoveParam != -1) trial->maxMoveParam = maxMoveParam;
    mc->initTrial(trial);
  } else if (typeStr.compare("pairMod") == 0) {
    shared_ptr<TrialPairMod> trial = make_shared<TrialPairMod>();
    if (maxMoveParam != -1) trial->maxMoveParam = maxMoveParam;
    mc->initTrial(trial);
  } else if (typeStr.compare("beta") == 0) {
    shared_ptr<TrialBeta> trial = make_shared<TrialBeta>();
    mc->initTrial(trial);
  } else if (typeStr.compare("xswap") == 0) {
    shared_ptr<TrialXSwap> trial = make_shared<TrialXSwap>();
    mc->initTrial(trial);

  // this is the non-legacy piece to keep when legacy is removed
  } else {
    shared_ptr<TrialTransform> trial = make_shared<TrialTransform>(type);
    if (maxMoveParam != -1) trial->maxMoveParam = maxMoveParam;
    mc->initTrial(trial);
  }
}
void transformTrial(shared_ptr<MC> mc, const char* type,
  double maxMoveParam = -1) {
  transformTrial(mc.get(), type, maxMoveParam);
}

/**
 * TrialCluster
 */
void clusterTrial(MC *mc, const char* type, const double clusterCut,
  const double maxMoveParam) {
  shared_ptr<TrialCluster> trial = make_shared<TrialCluster>(type);
  if (clusterCut != -1) trial->clusterCut = clusterCut;
  if (maxMoveParam != -1) trial->maxMoveParam = maxMoveParam;
  mc->initTrial(trial);
}
void clusterTrial(shared_ptr<MC> mc, const char* type, const double clusterCut,
  const double maxMoveParam) {
  clusterTrial(mc.get(), type, clusterCut, maxMoveParam);
}

void deleteTrial(MC *mc, const char* moltype) {
  shared_ptr<TrialDelete> trial = make_shared<TrialDelete>(moltype);
  mc->initTrial(trial);
}
void deleteTrial(shared_ptr<MC> mc, const char* moltype) {
  deleteTrial(mc.get(), moltype);
}
void deleteTrial(MC *mc) {
  shared_ptr<TrialDelete> trial = make_shared<TrialDelete>();
  mc->initTrial(trial);
}
void deleteTrial(shared_ptr<MC> mc) {
  deleteTrial(mc.get());
}

void addTrial(MC *mc, const char* moltype) {
  shared_ptr<TrialAdd> trial = make_shared<TrialAdd>(moltype);
  mc->initTrial(trial);
}
void addTrial(shared_ptr<MC> mc, const char* moltype) {
  addTrial(mc.get(), moltype);
}

void insertDeleteTrial(MC *mc, const char* moltype) {
  addTrial(mc, moltype);
  deleteTrial(mc, moltype);
}
void insertDeleteTrial(shared_ptr<MC> mc, const char* moltype) {
  insertDeleteTrial(mc.get(), moltype);
}

void gcaTrial(MC *mc, const int nMolTarg) {
  shared_ptr<TrialGCA> trial = make_shared<TrialGCA>();
  if (nMolTarg != -1) trial->targAcceptPer = nMolTarg;
  mc->initTrial(trial);
}
void gcaTrial(shared_ptr<MC> mc, const int nMolTarg) {
  gcaTrial(mc.get(), nMolTarg);
}

/**
 * PairMod
 */
void pairModTrial(MC *mc, const double maxMoveParam) {
  shared_ptr<TrialPairMod> trial = make_shared<TrialPairMod>();
  if (maxMoveParam != -1) trial->maxMoveParam = maxMoveParam;
  mc->initTrial(trial);
}
void pairModTrial(shared_ptr<MC> mc, const double maxMoveParam) {
  pairModTrial(mc.get(), maxMoveParam);
}
void pairModTrial(WLTMMC *mc, double maxMoveParam) {
  shared_ptr<TrialPairMod> trial = make_shared<TrialPairMod>();
  // by default, the max move parameter is the size of the WL bin
  if (maxMoveParam == -1) {
    maxMoveParam = mc->c()->mBin();
  }
  trial->maxMoveParam = maxMoveParam;
  mc->initTrial(trial);
}
void pairModTrial(shared_ptr<WLTMMC> mc, const double maxMoveParam) {
  pairModTrial(mc.get(), maxMoveParam);
}

/**
 * beta
 */
void betaTrial(MC *mc, const double maxMoveParam) {
  shared_ptr<TrialBeta> trial = make_shared<TrialBeta>();
  if (maxMoveParam != -1) trial->maxMoveParam = maxMoveParam;
  mc->initTrial(trial);
}
void betaTrial(shared_ptr<MC> mc, const double maxMoveParam) {
  betaTrial(mc.get(), maxMoveParam);
}
void betaTrial(WLTMMC *mc, double maxMoveParam) {
  shared_ptr<TrialBeta> trial = make_shared<TrialBeta>();
  // by default, the max move parameter is the size of the WL bin
  if (maxMoveParam == -1) {
    maxMoveParam = mc->c()->mBin();
  }
  trial->maxMoveParam = maxMoveParam;
  mc->initTrial(trial);
}
void betaTrial(shared_ptr<WLTMMC> mc, const double maxMoveParam) {
  betaTrial(mc.get(), maxMoveParam);
}

/**
 * grow
 */
void growTrial(MC *mc, const char* moltype, const int nStages) {
  shared_ptr<TrialGrow> trial = make_shared<TrialGrow>
    (mc->space(), mc->pair(), mc->criteria(), moltype, nStages);
  mc->initTrial(trial);
}
void growTrial(shared_ptr<MC> mc, const char* moltype, const int nStages) {
  growTrial(mc.get(), moltype, nStages);
}

/**
 * configbias
 */
void initConfigBias(MC *mc, const char* fileName, const int insDelFlag,
  const int dualCut) {
  shared_ptr<TrialConfigBias> trial = make_shared<TrialConfigBias>
    (mc->space(), mc->pair(), mc->criteria(), fileName);
  trial->insDelFlag = insDelFlag;
  trial->initDualCut(dualCut);
  trial->file2MC(fileName, mc);
}

void initConfigBias(shared_ptr<MC> mc, const char* fileName,
  const int insDelFlag, const int dualCut) {
  initConfigBias(mc.get(), fileName, insDelFlag, dualCut);
}

void xswapTrial(MC *mc) {
  shared_ptr<TrialXSwap> trial = make_shared<TrialXSwap>();
  mc->initTrial(trial);
}
void xswapTrial(shared_ptr<MC> mc) {
  xswapTrial(mc.get());
}

void swapTrial(MC *mc, const char* molTypeA, const char* molTypeB) {
  shared_ptr<TrialSwap> trial = make_shared<TrialSwap>(molTypeA, molTypeB);
  mc->initTrial(trial);
}
void swapTrial(shared_ptr<MC> mc, const char* molTypeA, const char* molTypeB) {
  swapTrial(mc.get(), molTypeA, molTypeB);
}

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_

#endif  // UI_ABBREVIATED_H_
