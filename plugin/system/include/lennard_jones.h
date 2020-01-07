
#ifndef FEASST_SYSTEM_LENNARD_JONES_H_
#define FEASST_SYSTEM_LENNARD_JONES_H_

#include <memory>
#include "system/include/model_two_body.h"
#include "math/include/constants.h"

namespace feasst {

/**
  The Lennard-Jones potential is given by

  \f$ U_{LJ} = 4\epsilon [ (\sigma/r)^{2\alpha} - (\sigma/r)^\alpha ] \f$,

  where \f$\alpha=6\f$ is assumed by this class for optimization.

 */
class LennardJones : public ModelTwoBody {
 public:
  LennardJones() {
    set_hard_sphere_threshold();
  }

  double energy(
      const double squared_distance,
      const int type1,
      const int type2,
      const ModelParams& model_params) const override {
    const double sigma = model_params.mixed_sigma()[type1][type2];
    const double sigma_squared = sigma*sigma;
    if (squared_distance < hard_sphere_threshold_sq_*sigma_squared) {
      return NEAR_INFINITY;
    }
    const double epsilon = model_params.mixed_epsilon()[type1][type2];
    const double rinv2 = sigma_squared/squared_distance;
    const double rinv6 = rinv2*rinv2*rinv2;
    const double en = 4.*epsilon*rinv6*(rinv6 - 1.);
    return en;
  }

  /// When the distance between sites does not exceed thresshold*sigma,
  /// then return NEAR_INFINITY energy.
  void set_hard_sphere_threshold(const double threshold = 0.2) {
    hard_sphere_threshold_sq_ = threshold*threshold; }

  /// Return the threshold for hard sphere interaction.
  double hard_sphere_threshold() const {
    return std::sqrt(hard_sphere_threshold_sq_); }
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
  const std::string class_name_ = "LennardJones";
  double hard_sphere_threshold_sq_;
};

inline std::shared_ptr<LennardJones> MakeLennardJones() {
  return std::make_shared<LennardJones>();
}

}  // namespace feasst

#endif  // FEASST_SYSTEM_LENNARD_JONES_H_