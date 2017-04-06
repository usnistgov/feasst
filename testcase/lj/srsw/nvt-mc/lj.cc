#include "pair_lj.h"
#include "mc.h"
#include "ui_abbreviated.h"

class AnalyzeMonkeyPatch : public Analyze {
 public:
  AnalyzeMonkeyPatch(Space *space, Pair *pair) : Analyze(space, pair) {}
  ~AnalyzeMonkeyPatch() {}
  Accumulator pe;
  void update() {
    pe.accumulate(pair_->peTot()/double(space_->nMol()));
  }
  void print() {
    cout << pe.average() << " +/- " << pe.blockStdev() << endl;
  }
};

int main() {
  Space space(3, 0);
  const double rho = 1e-3;
  const int nMol = 500;
  for (int dim = 0; dim < space.dimen(); ++dim) {
    space.lset(pow(double(nMol)/rho, 1./3.), dim);
  }
  stringstream molNameSS;
  molNameSS << space.install_dir() << "/forcefield/data.lj";
  space.addMolInit(molNameSS.str().c_str());
  PairLJ pair(&space, 3);   // potential truncation at 3
  pair.initEnergy();
  const double temperature = 0.9;
  CriteriaMetropolis criteria(1./temperature, 1.);
  MC mc(&space, &pair, &criteria);
  transformTrial(&mc, "translate", 0.1);
  mc.nMolSeek(nMol, molNameSS.str().c_str());
  mc.initLog("log", 1e4);
  mc.initMovie("movie", 1e4);
  AnalyzeMonkeyPatch an(&space, &pair);
  an.initFreq(1);
  an.initPrintFreq(1e5);
  mc.initAnalyze(&an);
  mc.runNumTrials(1e10);
}
