
#ifndef FEASST_MONTE_CARLO_TRIAL_SELECT_ANGLE_H_
#define FEASST_MONTE_CARLO_TRIAL_SELECT_ANGLE_H_

#include <memory>
#include "utils/include/arguments.h"
#include "monte_carlo/include/trial_select_bond.h"

namespace feasst {

/**
  A random particle of given type is selected if previously perturbed sites
    are not available.
  Select a single angle from two given anchor sites and one mobile site.
  The mobile site is directly bonded to the first anchor site, and the second
    anchor site is bonded to the first mobile site.
  The angle is defined as: anchor2 - anchor1 - mobile.
 */
class TrialSelectAngle : public TrialSelectBond {
 public:
  TrialSelectAngle(
    /**
      anchor_site2 : index of second anchor site.
     */
    const argtype& args = argtype()) : TrialSelectBond(args) {
    class_name_ = "TrialSelectAngle";
    Arguments args_(args);
    args_.dont_check();
    anchor_site2_ = args_.key("anchor_site2").integer();
  }

  // angle theta0 is added as a property
  // anchor is sized
  void precompute(System * system) override {
    TrialSelectBond::precompute(system);
    const Particle& part = system->configuration().particle_types().particle(particle_type());
    const int angle_type = part.angle(mobile_site(),
                                      anchor_site(),
                                      anchor_site2_).type();
    const Angle& angle = system->configuration().unique_types().particle(
      particle_type()).angle(angle_type);
    add_property("theta0", angle.property("theta0"));
    anchor_.add_site(0, anchor_site2_);
  }

  std::shared_ptr<TrialSelect> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit TrialSelectAngle(std::istream& istr);
  virtual ~TrialSelectAngle() {}

 protected:
  void serialize_trial_select_angle_(std::ostream& ostr) const;

 private:
  int anchor_site2_;
};

inline std::shared_ptr<TrialSelectAngle> MakeTrialSelectAngle(
    const argtype &args = argtype()) {
  return std::make_shared<TrialSelectAngle>(args);
}

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_TRIAL_SELECT_ANGLE_H_
