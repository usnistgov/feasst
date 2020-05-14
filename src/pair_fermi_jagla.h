/*
 * FEASST - Free Energy and Advanced Sampling Simulation Toolkit
 * http://pages.nist.gov/feasst, National Institute of Standards and Technology
 * Harold W. Hatch, harold.hatch@nist.gov
 *
 * Permission to use this data/software is contingent upon your acceptance of
 * the terms of LICENSE.txt and upon your providing
 * appropriate acknowledgments of NIST's creation of the data/software.
 */

#ifndef PAIR_FERMI_JAGLA_H_
#define PAIR_FERMI_JAGLA_H_

#include <memory>
#include <vector>
#include "./pair.h"

namespace feasst {

/**
 * Fermi-Jagla potential as described in:
 *  https://doi.org/10.1021/jp205098a
 *  https://doi.org/10.1039/C8SM00989A
 *
 * Note that, as described in https://doi.org/10.1039/C8SM00989A,
 * B_0 is the only parameter that is varied.
 * Thus, we set B_0 equal to Pair::epsij
 *
 * In addition, a hard sphere of diameter "r_s" avoids the singularity.
 *
 * Finally, by personal communication from Evan Pretti, some parameters
 * were changed to a slightly different level of percision than the previously
 * published values.
 */
class PairFermiJagla : public Pair {
 public:
  /// Constructor
  PairFermiJagla(Space* space, const argtype &args = argtype());

  // Potential parameters
  int n_c = 36;
  double epsilon_c = 10;
  double sigma_c = 0.2;
  double r_s = 0.8;
  double A_0 = 11.0346;   //11.035;
  double A_1 = 404.396;   //404.40;
  double A_2 = 1.0174094; //1.0174;
  double B_1 = 1044.5;
  double B_2 = 1.0305952; //1.0306;

  // Construct from restart file
  PairFermiJagla(Space* space, const char* fileName);
  virtual ~PairFermiJagla() {}
  virtual PairFermiJagla* clone(Space* space) const;

 protected:
  // See comments of derived class from Pair
  void pairSiteSite_(const int &iSiteType, const int &jSiteType,
    double * energy, double * force, int * neighbor, const double &dx,
    const double &dy, const double &dz);

  // defaults in constructor
  void defaultConstruction_();
};

/// Factory method
shared_ptr<PairFermiJagla> makePairFermiJagla(Space* space,
  const argtype &args = argtype());

}  // namespace feasst

#endif  // PAIR_FERMI_JAGLA_H_

