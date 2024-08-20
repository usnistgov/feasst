
#include "utils/test/utils.h"
#include "math/include/random_mt19937.h"
#include "math/include/histogram.h"
#include "configuration/include/group.h"
#include "system/include/hard_sphere.h"
#include "system/include/visit_model_cell.h"
#include "system/include/potential.h"
#include "system/include/thermo_params.h"
#include "monte_carlo/test/monte_carlo_utils.h"
#include "monte_carlo/include/monte_carlo.h"
#include "monte_carlo/include/trial_transfer.h"
#include "monte_carlo/include/trial_rotate.h"
#include "monte_carlo/include/trial_translate.h"
#include "monte_carlo/include/metropolis.h"
#include "monte_carlo/include/run.h"
#include "steppers/include/log.h"
#include "steppers/include/movie.h"
#include "steppers/include/check_energy.h"
#include "steppers/include/criteria_updater.h"
#include "steppers/include/criteria_writer.h"
#include "models/include/square_well.h"
#include "flat_histogram/include/transition_matrix.h"
#include "flat_histogram/include/wang_landau.h"
#include "flat_histogram/include/macrostate_num_particles.h"
#include "flat_histogram/include/flat_histogram.h"
#include "flat_histogram/include/clones.h"
#include "flat_histogram/include/window_exponential.h"
#include "patch/include/visit_model_inner_patch.h"
#include "patch/include/file_xyz_patch.h"
#include "patch/include/movie_patch.h"
#include "patch/include/solid_of_revolution_table.h"

namespace feasst {

std::unique_ptr<MonteCarlo> patchmc(const int min, const int max) {
  const double chi = 0.7;
  const double patch_angle_degrees = 2*std::asin(std::sqrt(chi/2))*180/PI;
  //DEBUG("patch_angle_degrees " << patch_angle_degrees);
  MonteCarlo mc;
  mc.set(MakeRandomMT19937({{"seed", "time"}}));
  { auto config = MakeConfiguration({{"cubic_side_length", "8"},
      {"patch_angle1", str(patch_angle_degrees)},
      {"particle_type0", install_dir() + "/plugin/patch/particle/two_patch_linear.fstprt"}});
      //{"particle_type0", install_dir() + "/plugin/patch/particle/janus.fstprt"}});
    config->add(MakeGroup({{"site_type0", "0"}}));
    mc.add(config);
  }
  mc.add(MakePotential(MakeHardSphere(),
                       MakeVisitModelCell({{"min_length", "1"}, {"cell_group_index", "1"}}),
                       {{"group_index", "1"}}));
  mc.add(MakePotential(
    MakeSquareWell(),
    MakeVisitModelCell(MakeVisitModelInnerPatch(),
        {{"min_length", "1.5"}, {"cell_group_index", "1"}}),
    {{"group_index", "1"}}));
  DEBUG(mc.configuration().model_params().select("patch_angle").str());
  DEBUG(mc.configuration().model_params().select("cos_patch_angle").str());
  DEBUG(mc.configuration().model_params().select("director").str());
  mc.set(MakeThermoParams({{"beta", str(1/0.7)}, {"chemical_potential", "-1.5"}}));
  //mc.set(MakeMetropolis());
  mc.set(MakeFlatHistogram(MakeMacrostateNumParticles({{"width", "1"}, {"max", feasst::str(max)}, {"min", feasst::str(min)}}),
    //MakeWangLandau({{"min_flatness", "100"}})));
    MakeTransitionMatrix({{"min_sweeps", "100"}})));
  mc.add(MakeTrialTranslate({{"tunable_param", "1"}}));
  mc.add(MakeTrialRotate({{"tunable_param", "40"}}));
  mc.add(MakeTrialTransfer({{"particle_type", "0"}, {"weight", "4"}}));
//  mc.add(MakeTrialAdd({{"particle_type", "0"}}));
//  mc.run(MakeRun({{"until_num_particles", "10"}}));
//  mc.run(MakeRemoveTrial({{"name", "TrialAdd"}}));
  const std::string trials_per = "1e5";
  mc.add(MakeLog({{"trials_per_write", trials_per}, {"output_file", "tmp/patch_nvt.txt"}}));
  mc.add(MakeMovie({{"trials_per_write", trials_per}, {"output_file", "tmp/patch_nvt.xyz"}}));
  mc.add(MakeCheckEnergy({{"trials_per_update", trials_per}}));
  mc.add(MakeMoviePatch({{"trials_per_write", trials_per}, {"output_file", "tmp/patch_nvt_vis.xyz"}}));
  mc.add(MakeCriteriaUpdater({{"trials_per_update", trials_per}}));
  mc.add(MakeCriteriaWriter({{"trials_per_write", trials_per}, {"output_file", "tmp/patch_fh.txt"}}));
  auto mc2 = test_serialize_unique(mc);
  //mc2.run_until_complete();
  return mc2;
}

TEST(MonteCarlo, patch) {
  auto mc2 = patchmc(0, 370);
  mc2->attempt(1e3);
  FileXYZPatch().write_for_vmd("tmp/test.xyz", mc2->configuration());
}

TEST(MonteCarlo, patch_LONG) {
  auto mc2 = patchmc(0, 370);
  mc2->attempt(1e6);
  FileXYZPatch().write_for_vmd("tmp/test.xyz", mc2->configuration());
}

//TEST(MonteCarlo, patch_clones_LONG) {
//  Clones clones;
//  clones.add(std::make_shared<MonteCarlo>(patchmc(0, 10)));
//  clones.add(std::make_shared<MonteCarlo>(patchmc(10, 20)));
//  DEBUG(clones.clone(0).configuration().model_params().select("patch_angle").value(1));
//  DEBUG(clones.clone(0).configuration().model_params().select("director").value(1));
//  DEBUG(clones.clone(0).configuration().model_params().select("cos_patch_angle").value(1));
//  clones.initialize(1);
////#  clones.get_clone(0)->initialize_criteria();
//}

//TEST(MonteCarlo, patch_clones2_LONG) {
//  std::vector<std::vector<int> > bounds = WindowExponential({
//    //{"maximum", "40"},
//    {"maximum", "20"},
//    //{"maximum", "370"},
//    {"minimum", "0"},
//    {"num", "2"},
//    //{"num", "4"},
//    {"extra_overlap", "0"},
//    {"alpha", "2"}}).boundaries();
//  Clones clones;
//  for (const std::vector<int> b : bounds) {
//    clones.add(std::make_shared<MonteCarlo>(patchmc(b[0], b[1])));
//  //clones.add(std::make_shared<MonteCarlo>(patchmc(10, 20)));
//  //clones.add(std::make_shared<MonteCarlo>(patchmc(10, 20)));
//  //clones.add(std::make_shared<MonteCarlo>(patchmc(10, 20)));
//  }
//  clones.initialize(1);
//  DEBUG("num " << clones.clone(0).configuration().num_particles());
//  DEBUG("num " << clones.clone(1).configuration().num_particles());
////  std::cout << "current 0 " << clones.clone(0).criteria().current_energy() << std::endl;
////  std::cout << "current 1 " << clones.clone(1).criteria().current_energy() << std::endl;
//  FileXYZPatch().write_for_vmd("tmp/test0.xyz", clones.clone(0).configuration());
//  FileXYZPatch().write_for_vmd("tmp/test1.xyz", clones.clone(1).configuration());
//  DEBUG("initialize 0");
//  //clones.get_clone(0)->initialize_criteria();
//  DEBUG("initialize 1");
//  //clones.get_clone(1)->initialize_criteria();
////  std::cout << "current 1 " << clones.clone(1).criteria().current_energy() << std::endl;
//  EXPECT_EQ(clones.clone(1).system().potential(1).group_index(), 1);
//  DEBUG(clones.clone(1).system().potential(1).model().class_name());
//  DEBUG(clones.clone(1).system().potential(1).visit_model().class_name());
////  std::cout << "current 0 " << clones.clone(0).criteria().current_energy() << std::endl;
////  std::cout << "current 1 " << clones.clone(1).criteria().current_energy() << std::endl;
//  //clones.initialize_and_run_until_complete();
//  clones.get_clone(0)->initialize_criteria();
////  FATAL("cosacut is changing for some reason, from the correct value of 0.3 to 0.866");
//}

TEST(MonteCarlo, patch_arglist) {
  auto mc = MakeMonteCarlo({{
    {"Configuration", {{"particle_type0", install_dir()+"/plugin/patch/particle/janus.fstprt"},
      {"xyz_file", "../plugin/patch/test/data/patch5.xyz"},
      {"cutoff", "3"},
      {"group0", "centers"}, {"centers_site_type0", "0"}}},
    {"Potential", {{"Model", "SquareWell"}, {"VisitModelInner", "VisitModelInnerPatch"}, {"group", "centers"}}},
    {"ThermoParams", {{"beta", "1"}}},
    {"Metropolis", {{}}}
  }});

  EXPECT_NEAR(-3., mc->criteria().current_energy(), NEAR_ZERO);
}

TEST(MonteCarlo, spherocylinder) {
  auto mc = MakeMonteCarlo({{
    {"Configuration", {{"particle_type0", install_dir()+"/plugin/patch/particle/spherocylinder.fstprt"},
      {"cubic_side_length", "8"}, {"group0", "centers"}, {"centers_site_type0", "0"}}},
    {"Potential", {{"Model", "SquareWell"}, {"VisitModelInner", "Spherocylinder"}, {"group", "centers"}}},
    {"ThermoParams", {{"beta", "1"}, {"chemical_potential0", "-1"}}},
    {"Metropolis", {{}}},
    {"TrialTranslate", {{"weight", "1."}, {"tunable_param", "1."}}},
    {"TrialRotate", {{"weight", "1."}, {"tunable_param", "1."}}},
    {"TrialAdd", {{"particle_type", "0"}}},
    {"Run", {{"until_num_particles", "20"}}},
    {"RemoveTrial", {{"name", "TrialAdd"}}},
    {"Checkpoint", {{"num_hours", "0.0001"}, {"checkpoint_file", "tmp/sphrst"}}},
    {"Log", {{"trials_per_write", str(1e0)}, {"output_file", "tmp/sph.txt"}}},
    {"Movie", {{"trials_per_write", str(1e0)}, {"output_file", "tmp/sph.xyz"}}},
    {"MovieSpherocylinder", {{"trials_per_write", str(1e0)}, {"output_file", "tmp/sphc.xyz"}}},
    {"CheckEnergy", {{"trials_per_update", str(1e0)}, {"tolerance", str(1e-9)}}},
    {"Tune", {{}}},
    {"Run", {{"num_trials", "1e2"}}},
  }});
}

TEST(MonteCarlo, SolidOfRevolution) {
  auto mc = MakeMonteCarlo({{
    //{"RandomMT19937", {{"seed", "123"}}},
    {"Configuration", {{"particle_type0", "../plugin/patch/particle/one_patch.fstprt"},
      {"cubic_side_length", "8"}, {"group0", "centers"}, {"centers_site_type0", "0"}, {"cutoff", "4"}}},
    {"Potential", {{"Model", "HardSphere"}, {"VisitModelInner", "SolidOfRevolutionTable"}, {"group", "centers"}, {"table_file", "../plugin/patch/test/data/tablek5l1.0d1.txt"}}},
    {"ThermoParams", {{"beta", "1"}, {"chemical_potential0", "-1"}}},
    {"Metropolis", {{}}},
    {"TrialTranslate", {{"weight", "1."}, {"tunable_param", "1."}}},
    {"TrialRotate", {{"weight", "1."}, {"tunable_param", "1."}}},
    {"TrialAdd", {{"particle_type", "0"}}},
    {"Run", {{"until_num_particles", "20"}}},
    {"RemoveTrial", {{"name", "TrialAdd"}}},
    {"Checkpoint", {{"num_hours", "0.0001"}, {"checkpoint_file", "tmp/sphrst"}}},
    {"Log", {{"trials_per_write", str(1e0)}, {"output_file", "tmp/sph.txt"}}},
    {"Movie", {{"trials_per_write", str(1e0)}, {"output_file", "tmp/sph.xyz"}}},
    {"CheckEnergy", {{"trials_per_update", str(1e0)}, {"tolerance", str(1e-9)}}},
    {"Tune", {{}}},
    {"Run", {{"num_trials", "1e2"}}},
  }});
}

}  // namespace feasst
