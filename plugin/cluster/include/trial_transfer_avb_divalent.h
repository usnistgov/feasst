
#ifndef FEASST_CLUSTER_TRIAL_ADD_AVB_DIVALENT_H_
#define FEASST_CLUSTER_TRIAL_ADD_AVB_DIVALENT_H_

#include <memory>
#include "utils/include/arguments.h"
#include "monte_carlo/include/trial.h"
#include "monte_carlo/include/trial_factory.h"

namespace feasst {

/**
  Attempt to add a particle with AVB as described in ComputeAddAVBDivalent.

  args:
  - particle_type_a: type of second added particle in AV of first.
  - site_index_a: index of site in type a that defines AV (default: 0).
  - particle_type_b: type of third added particle in AV of first.
  - site_index_b: index of site in type b that defines AV (default: 0).
 */
class TrialAddAVBDivalent : public Trial {
 public:
  TrialAddAVBDivalent(argtype args = argtype());
  TrialAddAVBDivalent(argtype * args);
  std::shared_ptr<Trial> create(std::istream& istr) const override {
    return std::make_shared<TrialAddAVBDivalent>(istr); }
  std::shared_ptr<Trial> create(argtype * args) const override {
    return std::make_shared<TrialAddAVBDivalent>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit TrialAddAVBDivalent(std::istream& istr);
  virtual ~TrialAddAVBDivalent() {}
};

inline std::shared_ptr<TrialAddAVBDivalent> MakeTrialAddAVBDivalent(argtype args = argtype()) {
  return std::make_shared<TrialAddAVBDivalent>(args); }

/**
  Attempt to add a particle with AVB as described in ComputeRemoveAVBDivalent.

  args:
  - particle_type_a: type of second added particle in AV of first.
  - site_index_a: index of site in type a that defines AV (default: 0).
  - particle_type_b: type of third added particle in AV of first.
  - site_index_b: index of site in type b that defines AV (default: 0).
 */
class TrialRemoveAVBDivalent : public Trial {
 public:
  TrialRemoveAVBDivalent(argtype args = argtype());
  TrialRemoveAVBDivalent(argtype * args);
  std::shared_ptr<Trial> create(std::istream& istr) const override {
    return std::make_shared<TrialRemoveAVBDivalent>(istr); }
  std::shared_ptr<Trial> create(argtype * args) const override {
    return std::make_shared<TrialRemoveAVBDivalent>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit TrialRemoveAVBDivalent(std::istream& istr);
  virtual ~TrialRemoveAVBDivalent() {}
};

inline std::shared_ptr<TrialRemoveAVBDivalent> MakeTrialRemoveAVBDivalent(argtype args = argtype()) {
  return std::make_shared<TrialRemoveAVBDivalent>(args); }

/// Attempt TrialAddAVBDivalent or TrialRemoveAVBDivalent with equal probability
class TrialTransferAVBDivalent : public TrialFactoryNamed {
 public:
  TrialTransferAVBDivalent(argtype args = argtype());
  TrialTransferAVBDivalent(argtype * args);
  std::shared_ptr<TrialFactoryNamed> create(argtype * args) const override {
    return std::make_shared<TrialTransferAVBDivalent>(args); }
  virtual ~TrialTransferAVBDivalent() {}
};

inline std::shared_ptr<TrialTransferAVBDivalent> MakeTrialTransferAVBDivalent(argtype args = argtype()) {
  return std::make_shared<TrialTransferAVBDivalent>(args); }

}  // namespace feasst

#endif  // FEASST_CLUSTER_TRIAL_ADD_AVB_DIVALENT_H_
