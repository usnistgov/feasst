
#ifndef FEASST_GROWTH_EXPANDED_TRIAL_MORPH_EXPANDED_H_
#define FEASST_GROWTH_EXPANDED_TRIAL_MORPH_EXPANDED_H_

#include <vector>
#include <string>
#include <memory>
#include "monte_carlo/include/trial.h"

namespace feasst {

class Random;

/**
  Grand canonical insert/deletion by gradual growth/shrink with expanded
  ensembles.

  As currently implemented, only one MacrostateMorph can be used at a time,
  due to the way that the current state is tracked.
 */
class TrialMorphExpanded : public Trial {
 public:
  /**
    Typically requires reference_index if multiple particles are to be morphed
    simultaneously.
   */
  TrialMorphExpanded(
    /// See MacrostateMorph for a description of the growth sequence.
    /// This growth sequence must be equivalent to MacrostateMorph.
    const std::vector<std::vector<int> > particle_type_growth_sequence,
    const argtype& args = argtype());
  void precompute(Criteria * criteria, System * system) override;
  bool attempt(Criteria * criteria, System * system, Random * random) override;
  std::shared_ptr<Trial> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit TrialMorphExpanded(std::istream& istr);
  virtual ~TrialMorphExpanded() {}
 private:
  std::vector<std::shared_ptr<Trial> > grow_;
  std::vector<std::shared_ptr<Trial> > shrink_;
  int current_state_;
};

inline std::shared_ptr<TrialMorphExpanded> MakeTrialMorphExpanded(
    const std::vector<std::vector<int> > particle_type_growth_sequence,
    const argtype &args = argtype()) {
  return std::make_shared<TrialMorphExpanded>(particle_type_growth_sequence,
                                              args);
}

}  // namespace feasst

#endif  // FEASST_GROWTH_EXPANDED_TRIAL_MORPH_EXPANDED_H_
