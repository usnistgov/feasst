#include "utils/test/utils.h"
#include "configuration/test/config_utils.h"
#include "system/include/hard_sphere.h"
#include "opt_lj/include/visit_model_opt_lj.h"

namespace feasst {

TEST(VisitModelOptLJ, reference_config) {
  Configuration config = lj_sample4();
  VisitModelOptLJ visit;
  visit.precompute(&config);
  HardSphere place_holder;
  Select one;
  one.add_particle(config.particle(0), 0);
  visit.compute(&place_holder, config.model_params(), one, &config);
  EXPECT_NEAR(-3.2639025245521616, visit.energy(), NEAR_ZERO);
}

}  // namespace feasst
