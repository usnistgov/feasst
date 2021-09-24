
#include "utils/test/utils.h"
#include "system/include/hard_sphere.h"
#include "monte_carlo/include/monte_carlo.h"
#include "monte_carlo/include/trials.h"
#include "monte_carlo/include/metropolis.h"
#include "monte_carlo/include/run.h"
#include "steppers/include/log_and_movie.h"
#include "models/include/square_well.h"
#include "patch/include/visit_model_inner_patch.h"
#include "patch/include/file_xyz_patch.h"
#include "patch/include/movie_patch.h"

namespace feasst {

TEST(MonteCarlo, patch_LONG) {
  MonteCarlo mc;
  { auto config = MakeConfiguration({{"cubic_box_length", "8"},
      {"particle_type0", install_dir() + "/plugin/patch/forcefield/two_patch_linear.fstprt"}});
      //{"particle_type0", install_dir() + "/plugin/patch/forcefield/janus.fstprt"}});
    config->add(MakeGroup({{"site_type0", "0"}}));
    mc.add(config);
  }
  mc.add(MakePotential(MakeSquareWell(),
                       MakeVisitModel(MakeVisitModelInnerPatch()),
                       {{"group_index", "1"}}));
  mc.add(MakePotential(MakeHardSphere(), {{"group_index", "1"}}));
  DEBUG(mc.configuration().model_params().select("patch_angle").str());
  DEBUG(mc.configuration().model_params().select("cos_patch_angle").str());
  DEBUG(mc.configuration().model_params().select("director").str());
  mc.set(MakeThermoParams({{"beta", "4"}, {"chemical_potential", "1"}}));
  mc.set(MakeMetropolis());
  mc.add(MakeTrialTranslate({{"tunable_param", "1"}}));
  mc.add(MakeTrialRotate({{"tunable_param", "40"}}));
  mc.add(MakeTrialAdd({{"particle_type", "0"}}));
  mc.perform(MakeRun({{"until_num_particles", "10"}}));
  mc.perform(MakeRemoveTrial({{"name", "TrialAdd"}}));
  const std::string steps_per = "1e4";
  mc.add(MakeLogAndMovie({{"steps_per", steps_per}, {"file_name", "tmp/patch_nvt"}}));
  mc.add(MakeMoviePatch({{"steps_per", steps_per}, {"file_name", "tmp/patch_nvt_vis.xyz"}}));
  mc.attempt(1e6);
  FileXYZPatch().write_for_vmd("tmp/test.xyz", mc.configuration());
}

}  // namespace feasst
