
#ifndef FEASST_CLUSTER_TRIAL_REMOVE_AVB_DIVALENT_H_
#define FEASST_CLUSTER_TRIAL_REMOVE_AVB_DIVALENT_H_

#include <memory>
#include "utils/include/arguments.h"
#include "monte_carlo/include/trial.h"
#include "monte_carlo/include/trial_factory.h"

namespace feasst {

/**
  Attempt to remove a particle of type "t" with site_index_t anywhere in Domain.
  Then, remove a second particle of type "a" with site_index_a in the AV of
  site_index_t.
  Finally, remove a third particle of type "a" with site_index_a in the AV of
  site_index_t.

  See TrialAddAVBDivalent for derivation of the acceptance probability that
  is the reverse of this Trial.

  There are modifications to make for this reverse move considering that the
  acceptance probability is computed before the removal takes place.

  The number of sites to select in the AV already contains the site added from
  the old state,

  \f$ N^{s,AV}_a + [1,2] \rightarrow N^{s,AV}_a \f$.

  The number of particles of type "t" already contains the first particle added
  from the old state.

  \f$ N_t + 1 \rightarrow N_t \f$
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
