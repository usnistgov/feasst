
#include "utils/test/utils.h"
#include "system/include/hard_sphere.h"
#include "system/include/visit_model_cell.h"
#include "monte_carlo/include/monte_carlo.h"
#include "monte_carlo/include/trials.h"
#include "monte_carlo/include/metropolis.h"
#include "monte_carlo/include/run.h"
#include "steppers/include/log_and_movie.h"
#include "steppers/include/check_energy.h"
#include "steppers/include/criteria_updater.h"
#include "steppers/include/criteria_writer.h"
#include "models/include/square_well.h"
#include "flat_histogram/include/transition_matrix.h"
#include "flat_histogram/include/wang_landau.h"
#include "flat_histogram/include/macrostate_num_particles.h"
#include "flat_histogram/include/flat_histogram.h"
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
  INFO("serializing after config");
  MonteCarlo mc3 = test_serialize(mc);
  mc.add(MakePotential(MakeHardSphere(),
                       MakeVisitModelCell({{"min_length", "1"}, {"cell_group", "1"}}),
                       {{"group_index", "1"}}));
  INFO("serializing after first potential");
  MonteCarlo mc4 = test_serialize(mc);
  mc.add(MakePotential(
    MakeSquareWell(),
    MakeVisitModelCell(MakeVisitModelInnerPatch(),
        {{"min_length", "1.5"}, {"cell_group", "1"}}),
    {{"group_index", "1"}}));
  INFO("serializing after second potential");
  MonteCarlo mc5 = test_serialize(mc);
  DEBUG(mc.configuration().model_params().select("patch_angle").str());
  DEBUG(mc.configuration().model_params().select("cos_patch_angle").str());
  DEBUG(mc.configuration().model_params().select("director").str());
  mc.set(MakeThermoParams({{"beta", str(1/0.7)}, {"chemical_potential", "-1.5"}}));
  //mc.set(MakeMetropolis());
  mc.set(MakeFlatHistogram(MakeMacrostateNumParticles({{"width", "1"}, {"max", "370"}, {"min", "0"}}),
    MakeWangLandau({{"min_flatness", "100"}})));
    //MakeTransitionMatrix({{"min_sweeps", "100"}})));
  mc.add(MakeTrialTranslate({{"tunable_param", "1"}}));
  mc.add(MakeTrialRotate({{"tunable_param", "40"}}));
  mc.add(MakeTrialTransfer({{"particle_type", "0"}, {"weight", "4"}}));
//  mc.add(MakeTrialAdd({{"particle_type", "0"}}));
//  mc.run(MakeRun({{"until_num_particles", "10"}}));
//  mc.run(MakeRemoveTrial({{"name", "TrialAdd"}}));
  INFO("serializing after trials");
  MonteCarlo mc6 = test_serialize(mc);
  const std::string steps_per = "1e6";
  mc.add(MakeLogAndMovie({{"steps_per", steps_per}, {"file_name", "tmp/patch_nvt"}}));
  mc.add(MakeCheckEnergy({{"steps_per", steps_per}}));
  INFO("serializing after check en");
  MonteCarlo m1 = test_serialize(mc);
  mc.add(MakeMoviePatch({{"steps_per", steps_per}, {"file_name", "tmp/patch_nvt_vis.xyz"}}));
  INFO("serializing after patch movie");
  MonteCarlo mc7 = test_serialize(mc);
  mc.add(MakeCriteriaUpdater({{"steps_per", steps_per}}));
  INFO("serializing after updater");
  MonteCarlo mc8 = test_serialize(mc);
  mc.add(MakeCriteriaWriter({{"steps_per", steps_per}, {"file_name", "tmp/patch_fh.txt"}}));
  INFO("serializing before run");
  MonteCarlo mc2 = test_serialize(mc);
  //mc2.run_until_complete();
  mc2.attempt(1e6);
  FileXYZPatch().write_for_vmd("tmp/test.xyz", mc2.configuration());
}

}  // namespace feasst
