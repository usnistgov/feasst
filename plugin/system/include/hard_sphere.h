
#ifndef FEASST_SYSTEM_HARD_SPHERE_H_
#define FEASST_SYSTEM_HARD_SPHERE_H_

#include <string>
#include <memory>
#include "system/include/model_two_body.h"

namespace feasst {

/**
  The potential energy for a hard sphere model of radius \f$\sigma\f$ and
  separation distance \f$r\f$ is as follows

  \f$ U_{HS} = \left\{
    \begin{array}{lr}
      \infty & r \le \sigma \\
      0 & r > \sigma
    \end{array}
  \right. \f$

  The cutoff model parameter must also be specified, or essentially equated to
  sigma in order to not skip this HardSphere interaction.
  For this reason, HardSphere-only models are often used in conjunction with the
  Configuration argument "set_cutoff_min_to_sigma".
 */
class HardSphere : public ModelTwoBody {
 public:
  HardSphere() { class_name_ = "HardSphere"; }

  double energy(
    const double squared_distance,
    const int type1,
    const int type2,
    const ModelParams& model_params) override;

  std::shared_ptr<Model> create(std::istream& istr) const override {
    return std::make_shared<HardSphere>(istr); }
  void serialize(std::ostream& ostr) const override;
  explicit HardSphere(std::istream& istr);
  virtual ~HardSphere() {}
};

inline std::shared_ptr<HardSphere> MakeHardSphere() {
  return std::make_shared<HardSphere>();
}

}  // namespace feasst

#endif  // FEASST_SYSTEM_HARD_SPHERE_H_
