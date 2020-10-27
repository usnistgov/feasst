
#ifndef FEASST_MONTE_CARLO_PERTURB_DISTANCE_H_
#define FEASST_MONTE_CARLO_PERTURB_DISTANCE_H_

#include "monte_carlo/include/perturb_move.h"

namespace feasst {

// HWH: could enable tuning and position w.r.t. previous bond placement
//      for higher acceptance probability
/**
 * Put first site in selection in a sphere about the first site in anchor.
 * For bond potentials, the equilibrium length and spring constant are as
 * described in Random::bond_length
 * Currently implemented for harmonic bonds (exponent: 2), but could add an
 * optional exponent model parameter to generalize this.
 */
class PerturbDistance : public PerturbMove {
 public:
  explicit PerturbDistance(const argtype& args = argtype());

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
};

inline std::shared_ptr<PerturbDistance> MakePerturbDistance(
    const argtype& args = argtype()) {
  return std::make_shared<PerturbDistance>(args);
}

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_PERTURB_DISTANCE_H_
