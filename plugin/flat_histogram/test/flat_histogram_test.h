#include "flat_histogram/include/flat_histogram.h"
#include "flat_histogram/include/macrostate_num_particles.h"
#include "flat_histogram/include/transition_matrix.h"
#include "flat_histogram/include/wang_landau.h"

namespace feasst {

inline std::shared_ptr<FlatHistogram> crit_fh(const int crit_type) {
  auto criteria = MakeFlatHistogram({
    {"beta", str(1./1.5)},
    {"chemical_potential", "-2.352321"}
  });
  criteria->set(MakeMacrostateNumParticles(
    Histogram({{"width", "1"}, {"max", "5"}}),
    {{"soft_max", "5"}}
  ));
  if (crit_type == 0) {
    criteria->set(MakeTransitionMatrix({
      {"min_sweeps", "10"},
      // {"num_steps_to_update", str(1e5)},  // benchmark 1.7 seconds
      {"num_steps_to_update", "1"}, // fast
    }));
  } else {
    criteria->set(MakeWangLandau({{"min_flatness", "1"}}));
  }
  return criteria;
}

}  // namespace feasst
