
#include "utils/test/utils.h"
#include "math/include/utils_math.h"
#include "configuration/include/configuration.h"
#include "configuration/include/group.h"
#include "system/include/hard_sphere.h"
#include "system/include/potential.h"
#include "system/include/thermo_params.h"
#include "system/include/visit_model.h"
#include "system/include/visit_model_inner.h"
#include "monte_carlo/include/monte_carlo.h"
#include "monte_carlo/include/trial_rotate.h"
#include "monte_carlo/include/trial_translate.h"
#include "models/include/square_well.h"
#include "mayer/include/mayer_sampling.h"
#include "patch/include/visit_model_inner_patch.h"

namespace feasst {

/// compare to Equation 8 of https://doi.org/10.1063/1.1569473
TEST(MayerSampling, patch_LONG) {
  for (const double degrees : {180, 90}) {
    MonteCarlo mc;
    const std::string fstprt = install_dir() + "/plugin/patch/particle/janus.fstprt";
    { auto config = MakeConfiguration({{"cubic_side_length", str(NEAR_INFINITY)}});
      config->add_particle_type(fstprt);
      config->add_particle_type(fstprt, "2");
      config->add(MakeGroup({{"site_type0", "0"}, {"site_type1", "2"}}));
      config->set_model_param("patch_angle", 1, degrees);
      config->set_model_param("patch_angle", 3, degrees);
      config->add_particle_of_type(0);
      config->add_particle_of_type(1);
      mc.add(config);
    }
    EXPECT_EQ(2, mc.configuration().num_particles());
    auto patch = MakeVisitModelInnerPatch();
    mc.add(MakePotential(MakeHardSphere(), {{"group_index", "1"}}));
    mc.add(MakePotential(MakeSquareWell(),
                         MakeVisitModel(patch),
                         {{"group_index", "1"}}));
    EXPECT_DOUBLE_EQ(patch->cos_patch_angle().value(0), 1);
    EXPECT_DOUBLE_EQ(patch->cos_patch_angle().value(1), std::cos(degrees_to_radians(degrees)));
    EXPECT_DOUBLE_EQ(patch->cos_patch_angle().value(2), 1);
    EXPECT_DOUBLE_EQ(patch->cos_patch_angle().value(3), std::cos(degrees_to_radians(degrees)));
    mc.add_to_reference(MakePotential(MakeHardSphere(), {{"group_index", "1"}}));
    const double beta = 0.1;
    mc.set(MakeThermoParams({{"beta", str(beta)}}));
    auto mayer = MakeMayerSampling();
    mc.set(mayer);
    INFO(mayer->current_energy());
    EXPECT_GT(mayer->current_energy(), 0);
    mc.add(MakeTrialTranslate({{"new_only", "true"}, {"reference_index", "0"}, {"tunable_param", "1"}, {"particle_type", "1"}}));
    mc.add(MakeTrialRotate({{"new_only", "true"}, {"reference_index", "0"}, {"tunable_param", "40"}}));
    const std::string trials_per = "1e4";
    //mc.add(MakeLog({{"trials_per_write", trials_per}, {"file_name", "tmp/patch.txt"}}));
    //mc.add(MakeMovie({{"trials_per_write", trials_per}, {"file_name", "tmp/patch.xyz"}}));
    mc.attempt(1e6);
    const double chi = std::pow(std::sin(degrees_to_radians(degrees)/2), 2);
    const double b2_reduced_analytical = 1-chi*chi*(std::pow(1.5, 3)-1)*(std::exp(beta)-1);
    INFO(mayer->mayer().str());
    INFO(mayer->mayer_ref().str());
    EXPECT_NEAR(b2_reduced_analytical, mayer->second_virial_ratio(), 10*mayer->second_virial_ratio_block_stdev());
  }
}

}  // namespace feasst
