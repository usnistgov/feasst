#ifndef FEASST_CHAIN_TRIAL_REPTATE_H_
#define FEASST_CHAIN_TRIAL_REPTATE_H_

#include <memory>
#include "utils/include/arguments.h"
#include "monte_carlo/include/trial_move.h"

namespace feasst {

/**
  Reptate a linear chain by taking one end and adding it to the other end.
  For heteropolymers, this perturbation changes the composition of all
  particles with the same type.
  Thus, individual heteropolymers should be added as unique particle types.
  The bond length is taken as the bond between site 0 and 1, and is assumed
  to be constant.
  Thus, as currently implemented, heteropolymers must have the same bond length.
  This trial may not be compatible with angle and dihedral potentials.
  Instead, use "reptate" in TrialGrow.
 */
class TrialReptate : public TrialMove {
 public:
  /**
    args:
    - Trial arguments.
    - SelectReptate arguments.
   */
  explicit TrialReptate(argtype args = argtype());
  explicit TrialReptate(argtype * args);
  std::shared_ptr<Trial> create(std::istream& istr) const override {
    return std::make_shared<TrialReptate>(istr); }
  std::shared_ptr<Trial> create(argtype * args) const override {
    return std::make_shared<TrialReptate>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit TrialReptate(std::istream& istr);
  virtual ~TrialReptate() {}
};

inline std::shared_ptr<TrialReptate> MakeTrialReptate(argtype args = argtype()) {
  return std::make_shared<TrialReptate>(args); }

}  // namespace feasst

#endif  // FEASST_CHAIN_TRIAL_REPTATE_H_
