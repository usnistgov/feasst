#include "flat_histogram/include/flat_histogram.h"
#include "flat_histogram/include/macrostate_num_particles.h"
#include "flat_histogram/include/wltm.h"

namespace feasst {

inline std::shared_ptr<FlatHistogram> crit_fh(
    /// 0 - transition matrix
    /// 1 - wang landau
    /// 2 - wltm
    const int crit_type) {
  auto criteria = MakeFlatHistogram({
    {"beta", str(1./1.5)},
    {"chemical_potential", "-2.352321"}
  });
  criteria->set(MakeMacrostateNumParticles(
    Histogram({{"width", "1"}, {"max", "5"}, {"min", "0"}})//,
//    {{"soft_max", "5"}, {"soft_min", "1"}}
  ));
  if (crit_type == 0) {
    criteria->set(MakeTransitionMatrix({
      {"min_sweeps", "10"},
//      // {"num_steps_to_update", str(1e5)},  // benchmark 1.7 seconds
//      {"num_steps_to_update", "1"}, // fast
    }));
  } else if (crit_type == 1) {
    criteria->set(MakeWangLandau({{"min_flatness", "20"}}));
  } else if (crit_type == 2) {
    criteria->set(MakeWLTM({
      {"collect_flatness", "15"},
      {"min_flatness", "20"},
      {"min_sweeps", "10"},
    }));
  } else {
    FATAL("unrecognized crit_type: " << crit_type);
  }
  return criteria;
}

}  // namespace feasst
