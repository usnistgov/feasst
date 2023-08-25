#include "utils/test/utils.h"
#include "monte_carlo/include/trial_select_particle.h"
#include "chain/include/perturb_library.h"

namespace feasst {

TEST(PerturbLibrary, serialize) {
  auto pert = MakePerturbLibrary({{"library_xyz",
    "../plugin/configuration/test/data/dimer4.xyz"}});
  PerturbLibrary pert2 = test_serialize(*pert);
}

TEST(PerturbLibrary, dimer) {
  System sys;
  sys.add(*MakeConfiguration({
    {"cubic_side_length", "8"},
    {"particle_type0", "../particle/dimer.fstprt"}}));
  sys.precompute();
  auto sel = MakeTrialSelectParticle({{"particle_type", "0"}});
  auto pert = MakePerturbLibrary({{"library_xyz",
    "../plugin/configuration/test/data/dimer4.xyz"}});
  sel->precompute(&sys);
  pert->precompute(sel.get(), &sys);
  PerturbLibrary pert2 = test_serialize(*pert);
  EXPECT_EQ(static_cast<int>(pert2.xyz().size()), 4);
  EXPECT_EQ(pert2.xyz()[0][0].coord(0), 0);
  EXPECT_EQ(pert2.xyz()[0][1].coord(0), 0);
  EXPECT_EQ(pert2.xyz()[1][0].coord(0), 0);
  EXPECT_EQ(pert2.xyz()[1][1].coord(0), 1);
  EXPECT_EQ(pert2.xyz()[2][0].coord(1), 0);
  EXPECT_EQ(pert2.xyz()[2][1].coord(1), 1);
  EXPECT_EQ(pert2.xyz()[3][0].coord(2), 0);
  EXPECT_EQ(pert2.xyz()[3][1].coord(2), 1);
}

}  // namespace feasst
