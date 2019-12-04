
#ifndef FEASST_MONTE_CARLO_PERTURB_DISTANCE_H_
#define FEASST_MONTE_CARLO_PERTURB_DISTANCE_H_

#include "monte_carlo/include/perturb_move.h"

namespace feasst {

// HWH: could enable tuning and position w.r.t. previous bond placement
//      for higher acceptance probability
/**
 * Put first site in selection in a sphere about the first site in anchor.
 */
class PerturbDistance : public PerturbMove {
 public:
  explicit PerturbDistance(const argtype& args = argtype());

  /// Compute and store the distance from the bond_length property in select.
  void precompute(TrialSelect * select, System * system) override;

  /// Return the distance.
  double distance() const { return distance_; }

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
};

inline std::shared_ptr<PerturbDistance> MakePerturbDistance(
    const argtype& args = argtype()) {
  return std::make_shared<PerturbDistance>(args);
}

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_PERTURB_DISTANCE_H_
