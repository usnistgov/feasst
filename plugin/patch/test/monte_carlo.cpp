
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
  const std::string trials_per = "1e5";
  auto mc = MakeMonteCarlo({{
    {"RandomMT19937", {{"seed", "123"}}},
    {"Configuration", {{"cubic_side_length", "8"},
      {"patch_angle1", str(patch_angle_degrees)},
      {"particle_type0", "../plugin/patch/particle/two_patch_linear.txt"},
      {"group0", "centers"}, {"centers_site_type0", "0"}}},
    {"Potential", {{"Model", "HardSphere"},
                   {"VisitModel", "VisitModelCell"}, {"min_length", "1"}, {"cell_group_index", "1"},
                   {"group_index", "1"}}},
    {"Potential", {{"Model", "SquareWell"}, {"VisitModel", "VisitModelCell"},
                   {"VisitModelInner", "VisitModelInnerPatch"},
                   {"min_length", "1.5"}, {"cell_group_index", "1"},
                   {"group_index", "1"}}},
    {"ThermoParams", {{"beta", str(1/0.7)}, {"chemical_potential", "-1.5"}}},
    {"FlatHistogram", {{"Macrostate", "MacrostateNumParticles"}, {"width", "1"}, {"max", feasst::str(max)}, {"min", feasst::str(min)},
                      {"Bias", "TransitionMatrix"}, {"min_sweeps", "100"}}},
    {"TrialTranslate", {{"tunable_param", "1"}}},
    {"TrialRotate", {{"tunable_param", "40"}}},
    {"TrialTransfer", {{"particle_type", "0"}, {"weight", "4"}}},
    {"Log", {{"trials_per_write", trials_per}, {"output_file", "tmp/patch_nvt.txt"}}},
    {"Movie", {{"trials_per_write", trials_per}, {"output_file", "tmp/patch_nvt.xyz"}}},
    {"CheckEnergy", {{"trials_per_update", trials_per}}},
    {"MoviePatch", {{"trials_per_write", trials_per}, {"output_file", "tmp/patch_nvt_vis.xyz"}}},
    {"CriteriaUpdater", {{"trials_per_update", trials_per}}},
    {"CriteriaWriter", {{"trials_per_write", trials_per}, {"output_file", "tmp/patch_fh.txt"}}},
  }}, true);
  auto mc2 = test_serialize_unique(*mc);
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

TEST(MonteCarlo, patch_arglist) {
  auto mc = MakeMonteCarlo({{
    {"Configuration", {{"particle_type0", "../plugin/patch/particle/janus.txt"},
      {"xyz_file", "../plugin/patch/test/data/patch5.xyz"},
      {"cutoff", "3"},
      {"group0", "centers"}, {"centers_site_type0", "0"}}},
    {"Potential", {{"Model", "SquareWell"}, {"VisitModelInner", "VisitModelInnerPatch"}, {"group", "centers"}}},
    {"ThermoParams", {{"beta", "1"}}},
    {"Metropolis", {{}}}
  }}, true);

  EXPECT_NEAR(-3., mc->criteria().current_energy(), NEAR_ZERO);
}

TEST(MonteCarlo, spherocylinder) {
  auto mc = MakeMonteCarlo({{
    {"Configuration", {{"particle_type0", "../plugin/patch/particle/spherocylinder.txt"},
      {"cubic_side_length", "8"}, {"group0", "centers"}, {"centers_site_type0", "0"}}},
    {"Potential", {{"Model", "SquareWell"}, {"VisitModelInner", "Spherocylinder"}, {"group", "centers"}}},
    {"ThermoParams", {{"beta", "1"}, {"chemical_potential0", "-1"}}},
    {"Metropolis", {{}}},
    {"TrialTranslate", {{"weight", "1."}, {"tunable_param", "1."}}},
    {"TrialRotate", {{"weight", "1."}, {"tunable_param", "1."}}},
    {"TrialAdd", {{"particle_type", "0"}}},
    {"Run", {{"until_num_particles", "20"}}},
    {"Remove", {{"name", "TrialAdd"}}},
    {"Checkpoint", {{"num_hours", "0.0001"}, {"checkpoint_file", "tmp/sphrst"}}},
    {"Log", {{"trials_per_write", str(1e0)}, {"output_file", "tmp/sph.txt"}}},
    {"Movie", {{"trials_per_write", str(1e0)}, {"output_file", "tmp/sph.xyz"}}},
    {"MovieSpherocylinder", {{"trials_per_write", str(1e0)}, {"output_file", "tmp/sphc.xyz"}}},
    {"CheckEnergy", {{"trials_per_update", str(1e0)}, {"tolerance", str(1e-9)}}},
    {"Tune", {{}}},
    {"Run", {{"num_trials", "1e2"}}},
  }}, true);
}

TEST(MonteCarlo, SolidOfRevolution) {
  auto mc = MakeMonteCarlo({{
    //{"RandomMT19937", {{"seed", "123"}}},
    {"Configuration", {{"particle_type0", "../plugin/patch/particle/one_patch.txt"},
      {"cubic_side_length", "8"}, {"group0", "centers"}, {"centers_site_type0", "0"}, {"cutoff", "4"}}},
    {"Potential", {{"Model", "HardSphere"}, {"VisitModelInner", "SolidOfRevolutionTable"}, {"group", "centers"}, {"table_file", "../plugin/patch/test/data/tablek5l1.0d1.txt"}}},
    {"ThermoParams", {{"beta", "1"}, {"chemical_potential0", "-1"}}},
    {"Metropolis", {{}}},
    {"TrialTranslate", {{"weight", "1."}, {"tunable_param", "1."}}},
    {"TrialRotate", {{"weight", "1."}, {"tunable_param", "1."}}},
    {"TrialAdd", {{"particle_type", "0"}}},
    {"Run", {{"until_num_particles", "20"}}},
    {"Remove", {{"name", "TrialAdd"}}},
    {"Checkpoint", {{"num_hours", "0.0001"}, {"checkpoint_file", "tmp/sphrst"}}},
    {"Log", {{"trials_per_write", str(1e0)}, {"output_file", "tmp/sph.txt"}}},
    {"Movie", {{"trials_per_write", str(1e0)}, {"output_file", "tmp/sph.xyz"}}},
    {"CheckEnergy", {{"trials_per_update", str(1e0)}, {"tolerance", str(1e-9)}}},
    {"Tune", {{}}},
    {"Run", {{"num_trials", "1e2"}}},
  }}, true);
}

}  // namespace feasst
