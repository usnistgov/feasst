
#ifndef FEASST_GROWTH_EXPANDED_TRIAL_MORPH_H_
#define FEASST_GROWTH_EXPANDED_TRIAL_MORPH_H_

#include <vector>
#include <string>
#include <memory>
#include "monte_carlo/include/trial.h"

namespace feasst {

/**
  Morph random particle(s) of given type into a different particle type(s).
  Typically requires the use of a reference index if multiple particles are to
  be morphed simultaneously.
 */
class TrialMorph : public Trial {
 public:
  /**
    args:
    - particle_type[i]: type of particle that will be morphed.
      The [i] is to be substituted for an integer 0, 1, 2, ...
    - particle_type_morph[i]: type of particle to morph into.
      The [i] is to be substituted for an integer 0, 1, 2, ...
      Each [i] should have a corresponding particle_type[i] argument.
   */
  TrialMorph(const argtype& args = argtype());
  std::shared_ptr<Trial> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit TrialMorph(std::istream& istr);
  virtual ~TrialMorph() {}
};

inline std::shared_ptr<TrialMorph> MakeTrialMorph(
    const argtype &args = argtype()) {
  return std::make_shared<TrialMorph>(args);
}

}  // namespace feasst

#endif  // FEASST_GROWTH_EXPANDED_TRIAL_MORPH_H_
