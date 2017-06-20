/**
 * \file
 *
 * \brief This is the user interface class
 *
 * The abbreviated user interface allows one-line instantiation of trials for
 * Monte Carlo. Note that while shared_ptr are encourged, one may also use raw
 * pointers to create the Trials; however, this could lead to memory leaks
 * and other issues. Shared pointers are not as simple to use for the python
 * interface.
 *
 */

#ifndef UI_ABBREVIATED_H_
#define UI_ABBREVIATED_H_

#include <memory>
#include "./mc.h"
#include "./mc_wltmmc.h"

#ifdef FEASST_NAMESPACE_
namespace feasst {
#endif  // FEASST_NAMESPACE_

void transformTrial(MC *mc, const char* type, double maxMoveParam = -1);
void transformTrial(shared_ptr<MC> mc, const char* type,
                    double maxMoveParam = -1);

void deleteTrial(MC *mc, const char* moltype);
void deleteTrial(shared_ptr<MC> mc, const char* moltype);
void deleteTrial(MC *mc);
void deleteTrial(shared_ptr<MC> mc);

void addTrial(MC *mc, const char* moltype);
void addTrial(shared_ptr<MC> mc, const char* moltype);

void growTrial(MC *mc, const char* moltype, const int nStages);
void growTrial(shared_ptr<MC> mc, const char* moltype, const int nStages);

void gcaTrial(MC *mc, const int nMolTarg = -1);
void gcaTrial(shared_ptr<MC> mc, const int nMolTarg = -1);

void pairModTrial(MC *mc, const double maxMoveParam = -1);
void pairModTrial(shared_ptr<MC> mc, const double maxMoveParam = -1);
void pairModTrial(WLTMMC *mc, double maxMoveParam = -1);
void pairModTrial(shared_ptr<WLTMMC> mc, const double maxMoveParam = -1);

void betaTrial(MC *mc, const double maxMoveParam = -1);
void betaTrial(shared_ptr<MC> mc, const double maxMoveParam = -1);
void betaTrial(WLTMMC *mc, double maxMoveParam = -1);
void betaTrial(shared_ptr<WLTMMC> mc, const double maxMoveParam = -1);

void clusterTrial(MC *mc, const char* type, const double clusterCut = -1,
                  const double maxMoveParam = -1);
void clusterTrial(shared_ptr<MC> mc, const char* type,
                  const double clusterCut = -1, const double maxMoveParam = -1);

void initConfigBias(MC *mc, const char* fileName, const int insDelFlag = 0,
  const int dualCut = 0);
void initConfigBias(shared_ptr<MC> mc, const char* fileName,
                    const int insDelFlag = 0, const int dualCut = 0);

void xswapTrial(MC *mc);
void xswapTrial(shared_ptr<MC> mc);

void swapTrial(MC *mc, const char* molTypeA, const char* molTypeB);
void swapTrial(shared_ptr<MC> mc, const char* molTypeA, const char* molTypeB);

#ifdef FEASST_NAMESPACE_
}  // namespace feasst
#endif  // FEASST_NAMESPACE_

#endif  // UI_ABBREVIATED_H_
