
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

  args:
  - particle_type[i]: type of particle that will be morphed.
    The [i] is to be substituted for an integer 0, 1, 2, ...
  - particle_type_morph[i]: type of particle to morph into.
    The [i] is to be substituted for an integer 0, 1, 2, ...
    Each [i] should have a corresponding particle_type[i] argument.
 */
std::shared_ptr<Trial> MakeTrialMorph(const argtype &args = argtype());

}  // namespace feasst

#endif  // FEASST_GROWTH_EXPANDED_TRIAL_MORPH_H_
