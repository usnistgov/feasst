
#ifndef FEASST_MONTE_CARLO_TRIAL_VOLUME_H_
#define FEASST_MONTE_CARLO_TRIAL_VOLUME_H_

#include <memory>
#include "utils/include/arguments.h"
#include "monte_carlo/include/trial.h"

namespace feasst {

/// Attempt to add a particle.
class TrialVolume : public Trial {
 public:
  explicit TrialVolume(const argtype& args = argtype());
  std::shared_ptr<Trial> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit TrialVolume(std::istream& istr);
  virtual ~TrialVolume() {}

 protected:
  void serialize_trial_volume_(std::ostream& ostr) const;
};

inline std::shared_ptr<TrialVolume> MakeTrialVolume(
    const argtype &args = argtype()) {
  return std::make_shared<TrialVolume>(args);
}

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_TRIAL_VOLUME_H_
