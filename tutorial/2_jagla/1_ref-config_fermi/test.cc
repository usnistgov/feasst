/**
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

#include "feasst.h"

int main() {  // JAGLA, REF_CONFIG_FERMI
  feasst::Space space(3);
  feasst::PairFermiJagla pair(&space, {{"rCut", "3"}});
  stringstream molNameA;
  molNameA << space.install_dir() << "/forcefield/data.lj";
  pair.initData(molNameA.str());
  stringstream molNameB;
  molNameB << space.install_dir() << "/forcefield/data.ljb";
  pair.initData(molNameB.str());

  // set the cut-off
  pair.equateRcutForAllTypes();

  // set the pairwise B_0 parameter
  pair.epsijset(0, 0, 1.);
  pair.epsijset(0, 1, 1.);
  pair.epsijset(1, 1, 1.);

  // create clones of Space and Pair to perform two separate tests
  feasst::Space * space2 = space.clone();
  feasst::PairFermiJagla * pair2 = pair.clone(space2);

  // first, test the interaction between two particles
  vector<double> xAdd(space.dimen(), 0.);
  pair.addMol(xAdd, molNameA.str().c_str());
  const double r = 1.0245;
  xAdd[0] = r;
  pair.addMol(xAdd, molNameB.str().c_str());
  pair.initEnergy();
  const double peExpected = -0.24869232467128419;
  //const double peExpected = -0.25082657931074182; // Fermi Jagla parameters changed
  ASSERT(fabs(pair.peTot() - peExpected) < feasst::DTOL,
    "analytical comparison pe: " << MAX_PRECISION << pair.peTot());

  // second, test a reference configuration
  pair2->epsijset(0, 0, 0.);
  pair2->epsijset(0, 1, 1.3218525572535642);
  pair2->epsijset(1, 1, 0.);
  std::stringstream config;
  config << space.install_dir() <<
    "/tutorial/2_jagla/1_ref-config_fermi/test_case.xyz";
  const int nA = 162;
  for (int i = 0; i < nA; ++i) pair2->addMol(molNameA.str());
  for (int i = 0; i < nA; ++i) pair2->addMol(molNameB.str());
  std::ifstream file(config.str().c_str());   // load the xyz file
  pair2->readXYZ(file);                       // read the xyz file
  pair2->initEnergy();                        // compute energy
  double pe = pair2->peTot()/space2->nMol();
  ASSERT(fabs(pe + 1.8207714765944178) < feasst::DTOL, "pe: " << MAX_PRECISION << pe);
  pair2->epsijset(0, 0, 1.);
  pair2->epsijset(1, 1, 1.);
  pair2->initEnergy();
  pe = pair2->peTot()/space2->nMol();
  ASSERT(fabs(pe + 2.0590968411085431) < feasst::DTOL, "pe: " << MAX_PRECISION << pe);
}
