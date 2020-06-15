
#ifndef FEASST_CLUSTER_TRIAL_ADD_AVB_DIVALENT_H_
#define FEASST_CLUSTER_TRIAL_ADD_AVB_DIVALENT_H_

#include <memory>
#include "utils/include/arguments.h"
#include "monte_carlo/include/trial.h"
#include "monte_carlo/include/trial_factory.h"
#include "system/include/neighbor_criteria.h"

namespace feasst {

/* HWH
for Add

for Del
stage0
  - TrialSelectParticle t
  - PerturbRemove
stage1
  - SelectParticleAVB a (from anchor)
  - PerturbRemoveAVB
stage2
  - same as a with b
 */

/// Attempt to add a particle with AVB as described in ComputeAddAVBDivalent.
class TrialAddAVBDivalent : public Trial {
 public:
  /**
    args:
    - particle_type_a: type of second added particle in AV of first.
    - site_index_a: index of site in type a that defines AV (default: 0).
    - particle_type_b: type of third added particle in AV of first.
    - site_index_b: index of site in type b that defines AV (default: 0).
   */
  TrialAddAVBDivalent(
    std::shared_ptr<NeighborCriteria> neighbor_criteria,
    const argtype& args = argtype());
  std::shared_ptr<Trial> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit TrialAddAVBDivalent(std::istream& istr);
  virtual ~TrialAddAVBDivalent() {}

 protected:
  void serialize_trial_add_avb_divalent_(std::ostream& ostr) const;
};

inline std::shared_ptr<TrialAddAVBDivalent> MakeTrialAddAVBDivalent(
    std::shared_ptr<NeighborCriteria> neighbor_criteria,
    const argtype &args = argtype()) {
  return std::make_shared<TrialAddAVBDivalent>(neighbor_criteria, args);
}

/// Attempt to add a particle with AVB as described in ComputeRemoveAVBDivalent.
class TrialRemoveAVBDivalent : public Trial {
 public:
  /**
    args:
    - particle_type_a: type of second added particle in AV of first.
    - site_index_a: index of site in type a that defines AV (default: 0).
    - particle_type_b: type of third added particle in AV of first.
    - site_index_b: index of site in type b that defines AV (default: 0).
   */
  TrialRemoveAVBDivalent(
    std::shared_ptr<NeighborCriteria> neighbor_criteria,
    const argtype& args = argtype());
  std::shared_ptr<Trial> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit TrialRemoveAVBDivalent(std::istream& istr);
  virtual ~TrialRemoveAVBDivalent() {}

 protected:
  void serialize_trial_remove_avb_divalent_(std::ostream& ostr) const;
};

inline std::shared_ptr<TrialRemoveAVBDivalent> MakeTrialRemoveAVBDivalent(
    std::shared_ptr<NeighborCriteria> neighbor_criteria,
    const argtype &args = argtype()) {
  return std::make_shared<TrialRemoveAVBDivalent>(neighbor_criteria, args);
}

/// Attempt TrialAddAVBDivalent or TrialRemoveAVBDivalent with equal probability
class TrialTransferAVBDivalent : public TrialFactory {
 public:
  explicit TrialTransferAVBDivalent(
    std::shared_ptr<NeighborCriteria> neighbor_criteria,
    const argtype& args = argtype());
};

inline std::shared_ptr<TrialTransferAVBDivalent> MakeTrialTransferAVBDivalent(
    std::shared_ptr<NeighborCriteria> neighbor_criteria,
    const argtype &args = argtype()) {
  return std::make_shared<TrialTransferAVBDivalent>(neighbor_criteria, args);
}

}  // namespace feasst

#endif  // FEASST_CLUSTER_TRIAL_ADD_AVB_DIVALENT_H_
