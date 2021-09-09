
#ifndef FEASST_MONTE_CARLO_PERTURB_DISTANCE_H_
#define FEASST_MONTE_CARLO_PERTURB_DISTANCE_H_

#include "system/include/rigid_bond.h"
#include "monte_carlo/include/perturb_move.h"

namespace feasst {

// HWH: could enable tuning and position w.r.t. previous bond placement
//      for higher acceptance probability
/**
  Put first site in selection in a spherical shell about the first site in
  anchor.
  For rigid bonds, the spherical shell has infinitesimal thickness (a sphere).
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

  /// Return the bond type.
  double bond_type() const { return bond_type_; }

  /// Return the randomly selected distance from the bond potential.
  double random_distance(const System& system,
    const TrialSelect* select,
    Random * random,
    double * bond_energy  // return the bond energy for Rosenbluth exclusion
  );

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
  int bond_type_ = -1;
  int potential_acceptance_;

  // temporary
  RigidBond bond_;

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
