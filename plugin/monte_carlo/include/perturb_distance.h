
#ifndef FEASST_MONTE_CARLO_PERTURB_DISTANCE_H_
#define FEASST_MONTE_CARLO_PERTURB_DISTANCE_H_

#include "monte_carlo/include/perturb_move.h"

namespace feasst {

// HWH: could enable tuning and position w.r.t. previous bond placement
//      for higher acceptance probability
/**
  Put first site in selection in a sphere about the first site in anchor.
  The following Bond Properties (e.g., /feasst/forcefield/README.rst):
  length, maximum_length and spring_constant and exponent are as described
  in Random::bond_length.
  The length in Bond Properties is equilibrium_length in Random::bond_length.
  If maximum_length is not provided, it is assumed to be twice the length.
  If spring_constant is not provided, the bond is assumed to be rigid.
  If exponent is not provided, it is assumed to be 2 (harmonic).
 */
class PerturbDistance : public PerturbMove {
 public:
  /**
    args:
    - potential_acceptance: index of intramolecular potential that will be used
      to select the move. Ignore if -1 (default: -1).
   */
  explicit PerturbDistance(argtype args = argtype());
  explicit PerturbDistance(argtype * args);

  /// Compute and store the distance from the bond_length property in select.
  /// Also store the spring constant.
  void precompute(TrialSelect * select, System * system) override;

  /// Return the equilibrium distance.
  double distance() const { return distance_; }

  /// Return the spring constant. If -1, the bond length is rigid.
  double spring_constant() const { return spring_constant_; }

  /// Return the randomly selected distance from the bond potential.
  double random_distance(Random * random,
    const double beta,  /// inverse temperature
    const int dimension) const;

  void move(System * system,
      TrialSelect * select,
      Random * random) override;

  // serialize
  std::shared_ptr<Perturb> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit PerturbDistance(std::istream& istr);
  virtual ~PerturbDistance() {}

 protected:
  void serialize_perturb_distance_(std::ostream& ostr) const;

 private:
  double distance_ = 1.;
  double spring_constant_ = -1;
  double maximum_length_ = -1;
  double exponent_ = -1;
  int potential_acceptance_;

  void move_once_(System * system,
      TrialSelect * select,
      Random * random);
};

inline std::shared_ptr<PerturbDistance> MakePerturbDistance(
    argtype args = argtype()) {
  return std::make_shared<PerturbDistance>(args);
}

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_PERTURB_DISTANCE_H_
