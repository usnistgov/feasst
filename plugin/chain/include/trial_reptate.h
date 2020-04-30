
#ifndef FEASST_CHAIN_TRIAL_REPTATE_H_
#define FEASST_CHAIN_TRIAL_REPTATE_H_

#include <string>
#include <memory>
#include "monte_carlo/include/trial_move.h"

namespace feasst {

/**
  Reptate a linear chain by taking one end and adding it to the other end.
  For heteropolymers, this perturbation changes the composition of all
  particles with the same type.
  Thus, individual heteropolymers should be added as unique particle types.
 */
class TrialReptate : public TrialMove {
 public:
  TrialReptate(
    /// These arguments are sent to both PerturbReptate and TrialStage.
    const argtype& args = argtype());
  std::shared_ptr<Trial> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit TrialReptate(std::istream& istr);
  virtual ~TrialReptate() {}
};

inline std::shared_ptr<TrialReptate> MakeTrialReptate(
    const argtype &args = argtype()) {
  return std::make_shared<TrialReptate>(args);
}

}  // namespace feasst

#endif  // FEASST_CHAIN_TRIAL_REPTATE_H_
