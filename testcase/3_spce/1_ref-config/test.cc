/**
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

#include "pair_lj_coul_ewald.h"

void expect_near(const double val1, const double val2, const double tolerance, const char* comment) {
  const double diff = val1 - val2;
  if (fabs(diff) > tolerance) {
    ASSERT(0, "Expected " << MAX_PRECISION << val1 << " to be the same as " << val2 << " within tolerance " << tolerance << " for case " << comment);
  }
}

int main() {  // SPCE, SRSW_REFCONF
  for (int srswConfig = 1; srswConfig <= 2; ++srswConfig) {
    feasst::Space space(3);
    space.lset(20.);
    feasst::PairLJCoulEwald pair(&space, 10.);
    stringstream molNameSS;
    molNameSS << space.install_dir() << "/forcefield/data.spce";
    pair.initData(molNameSS.str());
    pair.initKSpace(5.6,  // alpha*L
                    27);  // k^2 < k2max cutoff

    // Set atom-based cut-off for comparison to SRSW reference config
    // Warning: atom-based cut-offs will fail because an optimized routine for
    // molecule-based cutoff is used by Monte Carlo trials.
    pair.initAtomCut(1);
    for (int iType = 0; iType < space.nParticleTypes(); ++iType) {
      for (int jType = 0; jType < space.nParticleTypes(); ++jType) {
        pair.rCutijset(iType, jType, pair.rCut());
      }
    }

    // Read the configuration
    std::ostringstream ss;
    ss << space.install_dir() << "/testcase/3_spce/1_ref-config/spce_sample_config_periodic" << srswConfig << ".xyz";
    std::ifstream file(ss.str().c_str());
    pair.readxyz(file);
    pair.erfTable()->tableOff();  // turn off the error function table
    pair.initEnergy();

    // Compare with the SRSW
    const double KJperMol2Kelvin = 1e3/feasst::idealGasConstant;
    if (srswConfig == 1) {
      ASSERT(100 == space.nMol(), "number of molecules");
      expect_near(pair.peLJ()*KJperMol2Kelvin, 99538.736236886805, feasst::DTOL, "LJ");
      expect_near(pair.peLRC()*KJperMol2Kelvin, -823.71499511652178, feasst::DTOL, "LRC");
      expect_near(pair.peQReal()*KJperMol2Kelvin, -558888.919727851, 1e-9, "QReal");
      expect_near(pair.peQFrr()*KJperMol2Kelvin, 6270.0937839397402, 1e-11, "QFrr");
      // Note, this is the negative of the sum of E_self and E_intra in SRSW
      expect_near(pair.peQFrrSelf()*KJperMol2Kelvin, 34699.373693030466, 1e-9, "QFrrSelf");
      expect_near(pair.peTot()*KJperMol2Kelvin, -488603.17839517148, 1e-8, "Total");
    } else if (srswConfig == 2) {
      ASSERT(200 == space.nMol(), "number of molecules");
      expect_near(pair.peLJ()*KJperMol2Kelvin, 193712.42256683594, feasst::DTOL, "LJ");
      expect_near(pair.peLRC()*KJperMol2Kelvin, -3294.8599804660735, feasst::DTOL, "LRC");
      expect_near(pair.peQReal()*KJperMol2Kelvin, -1192948.6940245957, 2e-9, "QReal");
      expect_near(pair.peQFrr()*KJperMol2Kelvin, 6034.9503269097158, 1e-11, "QFrr");
      // Note, this is the negative of the sum of E_self and E_intra in SRSW
      expect_near(pair.peQFrrSelf()*KJperMol2Kelvin, 69398.747386013289, 1e-9, "QFrrSelf");
      expect_near(pair.peTot()*KJperMol2Kelvin, -1065894.9284973294, 2e-9, "Total");
    }
  }
}
