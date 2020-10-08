
#ifndef FEASST_SYSTEM_LENNARD_JONES_H_
#define FEASST_SYSTEM_LENNARD_JONES_H_

#include <string>
#include <memory>
#include "system/include/model_two_body.h"

namespace feasst {

/**
  The Lennard-Jones potential is given by

  \f$ U_{LJ} = 4\epsilon [ (\sigma/r)^{2\alpha} - (\sigma/r)^\alpha ] \f$,

  where \f$\alpha=6\f$ is assumed by this class for optimization.

 */
class LennardJones : public ModelTwoBody {
 public:
  LennardJones();

  double energy(
      const double squared_distance,
      const int type1,
      const int type2,
      const ModelParams& model_params) const override;

  /// When the distance between sites does not exceed thresshold*sigma,
  /// then return NEAR_INFINITY energy.
  void set_hard_sphere_threshold(const double threshold = 0.2) {
    hard_sphere_threshold_sq_ = threshold*threshold; }

  /// Return the threshold for hard sphere interaction.
  double hard_sphere_threshold() const;
  const double& hard_sphere_threshold_sq() const {
    return hard_sphere_threshold_sq_; }

  std::shared_ptr<Model> create(std::istream& istr) const override {
    return std::make_shared<LennardJones>(istr); }
  void serialize(std::ostream& ostr) const override;
  explicit LennardJones(std::istream& istr);
  virtual ~LennardJones() {}

 protected:
  void serialize_lennard_jones_(std::ostream& ostr) const;

 private:
  double hard_sphere_threshold_sq_;
};

inline std::shared_ptr<LennardJones> MakeLennardJones() {
  return std::make_shared<LennardJones>();
}

}  // namespace feasst

#endif  // FEASST_SYSTEM_LENNARD_JONES_H_
