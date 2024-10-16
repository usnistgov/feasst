
#ifndef FEASST_CHAIN_GHOST_TRIAL_GROW_H_
#define FEASST_CHAIN_GHOST_TRIAL_GROW_H_

#include "math/include/accumulator.h"
#include "math/include/histogram.h"
#include "monte_carlo/include/modify.h"
#include "monte_carlo/include/trial.h"
#include "monte_carlo/include/trial_factory.h"
#include "monte_carlo/include/always_reject.h"

namespace feasst {

/**
  Compute the Metropolis acceptance for a TrialGrow that is never accepted.

  This trial may not be currently usable to obtain Widom insertions of flexible
  molecules due to the exclusion of bonded energy terms.

  This class effectively functions as an Analyze, because it does not change the
  System.
  The Trial is always rejected; however, because the Trial interface is used
  then the System is temporarily modified.
 */
class GhostTrialGrow : public Modify {
 public:
  GhostTrialGrow() : Modify() {} // only use for deserialize_map.

  //@{
  /** @name Arguments
    - grow_file: file name for use with TrialGrowFile.
    - Stepper arguments.
   */
  explicit GhostTrialGrow(argtype args);
  explicit GhostTrialGrow(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  std::string header(const Criteria& criteria,
    const System& system,
    const TrialFactory& trials) const override;

  void initialize(Criteria * criteria,
    System * system,
    TrialFactory * trial_factory) override;

  void update(Criteria * criteria,
    System * system,
    Random * random,
    TrialFactory * trial_factory) override;

  std::string write(Criteria * criteria,
    System * system,
    TrialFactory * trial_factory) override;

  // serialize
  std::string class_name() const override { return std::string("GhostTrialGrow"); }
  std::shared_ptr<Modify> create(std::istream& istr) const override {
    return std::make_shared<GhostTrialGrow>(istr); }
  std::shared_ptr<Modify> create(argtype * args) const override {
    return std::make_shared<GhostTrialGrow>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit GhostTrialGrow(std::istream& istr);

  //@}
 private:
  TrialFactory grow_;
  AlwaysReject criteria_;
  Accumulator metropolis_prob_;
};

inline std::shared_ptr<GhostTrialGrow> MakeGhostTrialGrow(
    argtype args) {
  return std::make_shared<GhostTrialGrow>(args);
}

}  // namespace feasst

#endif  // FEASST_CHAIN_GHOST_TRIAL_GROW_H_
