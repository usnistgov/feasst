
#ifndef FEASST_CLUSTER_TRIAL_REMOVE_AVB_DIVALENT_H_
#define FEASST_CLUSTER_TRIAL_REMOVE_AVB_DIVALENT_H_

#include <memory>
#include "utils/include/arguments.h"
#include "monte_carlo/include/trial.h"
#include "monte_carlo/include/trial_factory.h"

namespace feasst {

/**
  Attempt to add a particle with AVB as described in ComputeRemoveAVBDivalent.
 */
class TrialRemoveAVBDivalent : public Trial {
 public:
  //@{
  /** @name Arguments
    - particle_type_a: type of second added particle in AV of first.
    - site_index_a: index of site in type a that defines AV (default: 0).
    - particle_type_b: type of third added particle in AV of first.
    - site_index_b: index of site in type b that defines AV (default: 0).
    - SelectParticleAVBDivalent arguments.
   */
  explicit TrialRemoveAVBDivalent(argtype args = argtype());
  explicit TrialRemoveAVBDivalent(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  std::shared_ptr<Trial> create(std::istream& istr) const override {
    return std::make_shared<TrialRemoveAVBDivalent>(istr); }
  std::shared_ptr<Trial> create(argtype * args) const override {
    return std::make_shared<TrialRemoveAVBDivalent>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit TrialRemoveAVBDivalent(std::istream& istr);
  virtual ~TrialRemoveAVBDivalent() {}
  //@}
};

inline std::shared_ptr<TrialRemoveAVBDivalent> MakeTrialRemoveAVBDivalent(argtype args = argtype()) {
  return std::make_shared<TrialRemoveAVBDivalent>(args); }

}  // namespace feasst

#endif  // FEASST_CLUSTER_TRIAL_REMOVE_AVB_DIVALENT_H_
