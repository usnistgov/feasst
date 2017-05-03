#include "pair_tabular.h"
#include <string>

namespace feasst {

PairTabular::PairTabular(Space* space) : Pair(space, 0.) {
}
PairTabular::PairTabular(Space* space,
         const char* fileName)
  : Pair(space, fileName) {
}

/**
 * initialize tables
 */
void PairTabular::initTableHard(const int iType, const int jType,
                                const char* fileName) {
  tabHardFileName_[iType][jType] = fileName;
  tabHard_[iType][jType] = make_shared<Table>(fileName);
  tabHard_[jType][iType] = tabHard_[iType][jType];
}
void PairTabular::initTableCut(const int iType, const int jType,
                               const char* fileName) {
  tabCutFileName_[iType][jType] = fileName;
  tabCut_[iType][jType] = make_shared<Table>(fileName);
  tabCut_[jType][iType] = tabCut_[iType][jType];
}
void PairTabular::initTablePE(const int iType, const int jType,
                              const char* fileName) {
  tabPEFileName_[iType][jType] = fileName;
  tabPE_[iType][jType] = make_shared<Table>(fileName);
  tabPE_[jType][iType] = tabPE_[iType][jType];
}

/**
 * initialize whether to include more than hard interactions
 */
void PairTabular::initHard(const int iType, const int jType,
                           const int flag) {
  if (flag == 0) {
    hardFlag_[iType][jType] = 0;
    hardFlag_[jType][iType] = 0;
  } else {
    hardFlag_[iType][jType] = 1;
    hardFlag_[jType][iType] = 1;
  }
}

}  // namespace feasst



