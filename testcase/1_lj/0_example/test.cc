#include "feasst.h"

int main() {  // LJ, NVT_EXAMPLE
  auto space = feasst::makeSpace(
    {{"dimen", "3"},              // 3D space
     {"boxLength", "8"}});        // cubic periodic boundaries
  auto pair = feasst::makePairLJ(space,
    {{"rCut", "3"},               // potential truncation
     {"cutType", "lrc"}});        // long range corrections
  auto criteria = feasst::makeCriteriaMetropolis(
    {{"beta", "1.2"}});           // beta = 1/k_B/T
  feasst::MC mc(pair, criteria);
  feasst::addTrialTransform(&mc,
    {{"transType", "translate"},       // attempt particle translations
     {"maxMoveParam", "0.1"}});   // maximum displacement for each dimension
  mc.nMolSeek(50);                // add particles
  mc.initLog("log", 1e4);         // output instantaneous values
  mc.initMovie("movie", 1e4);     // output xyz trajectory
  mc.runNumTrials(1e6);           // perform MC trials
}
