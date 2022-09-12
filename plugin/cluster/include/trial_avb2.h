
#ifndef FEASST_CLUSTER_TRIAL_AVB2_H_
#define FEASST_CLUSTER_TRIAL_AVB2_H_

#include <memory>
#include "utils/include/arguments.h"
#include "monte_carlo/include/trial_factory.h"

namespace feasst {

/**
  Attempt AVB2 in->out or out->in
  Only implemented for single-site particles
  TrialAVB2 below is recommended in most cases to ensure detailed-balance is
  satisfied.
  But this is used in special cases like Prefetch when avoiding TrialFactory.

  args:
    - out_to_in: if true, use out->in. Otherwise, in->out.
    - neighbor_index: NeighborCriteria index contained in System (default: 0).
 */
class TrialAVB2Half : public Trial {
 public:
  explicit TrialAVB2Half(argtype args = argtype());
  explicit TrialAVB2Half(argtype * args);
  std::shared_ptr<Trial> create(std::istream& istr) const override {
    return std::make_shared<TrialAVB2Half>(istr); }
  std::shared_ptr<Trial> create(argtype * args) const override {
    return std::make_shared<TrialAVB2Half>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit TrialAVB2Half(std::istream& istr);
  virtual ~TrialAVB2Half() {}
};

inline std::shared_ptr<TrialAVB2Half> MakeTrialAVB2Half(argtype args = argtype()) {
  return std::make_shared<TrialAVB2Half>(args); }

/// Attempt AVB2 in->out and out->in with equal probability.
/// Only implemented for single-site particles
/// See ComputeAVB2 for more information.
class TrialAVB2 : public TrialFactoryNamed {
 public:
  explicit TrialAVB2(argtype args = argtype());
  explicit TrialAVB2(argtype * args);
  std::shared_ptr<TrialFactoryNamed> create(argtype * args) const override {
    return std::make_shared<TrialAVB2>(args); }
  virtual ~TrialAVB2() {}
};

inline std::shared_ptr<TrialAVB2> MakeTrialAVB2(argtype args = argtype()) {
  return std::make_shared<TrialAVB2>(args); }

// Process AVB2 args, which can also be used in TrialGrow
void gen_avb2_args_(const bool out_to_in, argtype * args, argtype * perturb_args);

}  // namespace feasst

#endif  // FEASST_CLUSTER_TRIAL_AVB2_H_
