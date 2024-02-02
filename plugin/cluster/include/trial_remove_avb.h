
#ifndef FEASST_CLUSTER_TRIAL_REMOVE_AVB_H_
#define FEASST_CLUSTER_TRIAL_REMOVE_AVB_H_

#include <memory>
#include "utils/include/arguments.h"
#include "monte_carlo/include/trial_factory.h"

namespace feasst {

/**
  Attempt to remove a particle of type "a" with a bias for site_index_a
  to be in the aggrevation volume (AV) of site_index_t ("t" for target) of
  particle of type "t".

  See TrialAddAVB for derivation of the acceptance probability that is
  the reverse of this Trial.

  There are modifications to make for this reverse move considering that the
  acceptance probability is computed before the removal takes place.

  The number of sites to select in the AV already contains the site added from
  the old state,

  \f$ N^{s,AV}_a + 1 \rightarrow N^{s,AV}_a \f$.
 */
class TrialRemoveAVB : public Trial {
 public:
  explicit TrialRemoveAVB(argtype args = argtype());
  explicit TrialRemoveAVB(argtype * args);
  std::shared_ptr<Trial> create(std::istream& istr) const override {
    return std::make_shared<TrialRemoveAVB>(istr); }
  std::shared_ptr<Trial> create(argtype * args) const override {
    return std::make_shared<TrialRemoveAVB>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit TrialRemoveAVB(std::istream& istr);
  virtual ~TrialRemoveAVB() {}
};

inline std::shared_ptr<TrialRemoveAVB> MakeTrialRemoveAVB(argtype args = argtype()) {
  return std::make_shared<TrialRemoveAVB>(args); }

}  // namespace feasst

#endif  // FEASST_CLUSTER_TRIAL_REMOVE_AVB_H_
