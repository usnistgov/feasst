#include "utils/test/utils.h"
#include "monte_carlo/include/trial_select_particle.h"

namespace feasst {

TEST(TrialSelectParticle, serialize) {
  TrialSelectParticle add;
//  std::stringstream ss;
//  add.serialize(ss);
//  INFO(ss.str());
//  TrialSelectParticle add2(ss);
//  add2.serialize(ss);
//  INFO(ss.str());
  TrialSelectParticle add2 = test_serialize(add);
}

}  // namespace feasst
