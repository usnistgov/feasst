
#ifndef FEASST_GIBBS_CHECK_CONSTANT_VOLUME_H_
#define FEASST_GIBBS_CHECK_CONSTANT_VOLUME_H_

#include <map>
#include <string>
#include "monte_carlo/include/modify.h"

namespace feasst {

typedef std::map<std::string, std::string> argtype;

/**
  Check that the total volume of all Configurations, stored from the last
  update, is within the tolerance of the current total volume.
 */
class CheckConstantVolume : public ModifyUpdateOnly {
 public:
  //@{
  /** @name Arguments
    - tolerance: relative absolute difference between current total volume
      and last recorded total volume (default: 1e-4).
    - Stepper arguments.
  */
  explicit CheckConstantVolume(argtype args = argtype());
  explicit CheckConstantVolume(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  void update(Criteria * criteria,
    System * system,
    Random * random,
    TrialFactory * trial_factory) override;

  std::string class_name() const override { return std::string("CheckConstantVolume"); }

  // serialize
  void serialize(std::ostream& ostr) const override;
  std::shared_ptr<Modify> create(std::istream& istr) const override {
    return std::make_shared<CheckConstantVolume>(istr); }
  std::shared_ptr<Modify> create(argtype * args) const override {
    return std::make_shared<CheckConstantVolume>(args); }
  explicit CheckConstantVolume(std::istream& istr);

  //@}
 private:
  double tolerance_;
  double last_total_volume_ = -1.;
};

inline std::shared_ptr<CheckConstantVolume> MakeCheckConstantVolume(
    argtype args = argtype()) {
  return std::make_shared<CheckConstantVolume>(args);
}

}  // namespace feasst

#endif  // FEASST_GIBBS_CHECK_CONSTANT_VOLUME_H_
