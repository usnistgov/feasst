#include "utils/test/utils.h"
#include "math/include/utils_math.h"
#include "configuration/include/group.h"
#include "system/include/hard_sphere.h"
#include "system/include/potential.h"
#include "system/include/visit_model.h"
#include "system/include/visit_model_inner.h"
#include "monte_carlo/include/monte_carlo.h"
#include "models/include/square_well.h"
#include "mayer/include/mayer_sampling.h"
#include "patch/include/visit_model_inner_patch.h"
#include "patch/include/file_xyz_patch.h"

namespace feasst {

TEST(FileXYZPatch, serialize) {
  auto patch = std::make_shared<FileXYZPatch>(argtype({{"append", "true"}}));
  FileXYZPatch patch2 = test_serialize(*patch);
  EXPECT_TRUE(patch2.append());
}

TEST(FileXYZPatch, patch) {
  MonteCarlo mc;
  const std::string fstprt = "../plugin/patch/particle/janus.txt";
  { auto config = MakeConfiguration({{"cubic_side_length", "8"}});
    config->add_particle_type(fstprt);
    config->add_particle_type(fstprt);
    config->add(MakeGroup({{"site_type", "0,2"}}));
    config->add_particle_of_type(0);
    config->add_particle_of_type(1);
    mc.add(config);
  }
  EXPECT_EQ(2, mc.configuration().num_particles());
  mc.add(MakePotential(MakeSquareWell(),
                       MakeVisitModel(std::make_shared<VisitModelInnerPatch>()),
                       {{"group_index", "1"}}));
  FileXYZPatch().write_for_vmd("tmp/test.xyz", mc.configuration());
}

}  // namespace feasst
