#include "utils/test/utils.h"
#include "monte_carlo/include/trial_select_dihedral.h"

namespace feasst {

TEST(TrialSelectDihedral, serialize) {
  TrialSelectDihedral add({{"mobile_site", "3"}, {"anchor_site", "2"}, {"anchor_site2", "1"}, {"anchor_site3", "0"}});
//  std::stringstream ss;
//  add.serialize(ss);
//  INFO(ss.str());
//  TrialSelectDihedral add2(ss);
//  add2.serialize(ss);
//  INFO(ss.str());
  TrialSelectDihedral add3 = test_serialize(add);
}

}  // namespace feasst
