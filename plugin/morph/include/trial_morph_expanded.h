
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

  While this methodology is essentially equivalent to the growth expanded
  ensemble, it was renamed morph because it is more general than adding sites
  to a polymer chain.
  However, each particle type could be a polymer chain with a variable number
  of interacting sites (e.g., modified epsilon/sigma).
  Thus morph expanded is equivalent to growth expanded in this special case.

  Each particle type in the growth sequence can have simultaneous changes in
  ModelParams such as Epsilon, Sigma, CutOff and Charge as the particle
  is slowly morphed from an particle with a high insertion/deletion acceptance
  to the desired particle.
  Note that the chemical potential of each particle type may be manually-tuned
  in an attempt to reduce the free-energy differences during morph trials.
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
    argtype args = argtype());
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
    argtype args = argtype()) {
  return std::make_shared<TrialMorphExpanded>(particle_type_growth_sequence,
                                              args);
}

}  // namespace feasst

#endif  // FEASST_GROWTH_EXPANDED_TRIAL_MORPH_EXPANDED_H_
