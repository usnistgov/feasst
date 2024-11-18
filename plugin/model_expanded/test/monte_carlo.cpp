
#include "utils/test/utils.h"
#include "monte_carlo/include/monte_carlo.h"

namespace feasst {

TEST(MonteCarlo, TrialModel) {
  auto mc = MakeMonteCarlo({{
    //{"RandomMT19937", {{"seed", "123"}}},
    {"Configuration", {{"particle_type0", "../particle/lj.fstprt"},
      {"cubic_side_length", "8"}}},
    //{"Potential", {{"Model", "LennardJones"}}},
    {"Potential", {{"Model", "ModelExpanded"}, {"model_file", "../plugin/model_expanded/test/data/models.txt"}}},
    {"ThermoParams", {{"beta", "1"}, {"chemical_potential0", "-1"}}},
    {"Metropolis", {{}}},
    {"TrialTranslate", {{"weight", "1."}, {"tunable_param", "1."}}},
    {"TrialAdd", {{"particle_type", "0"}}},
    {"Run", {{"until_num_particles", "20"}}},
    {"RemoveTrial", {{"name", "TrialAdd"}}},
    {"Checkpoint", {{"num_hours", "0.0001"}, {"checkpoint_file", "tmp/merst"}}},
    //{"Metropolis", {{}}},
    {"FlatHistogram", {{"Macrostate", "MacrostateModel"}, {"width", "1"}, {"max", "1"}, {"Bias", "WangLandau"}, {"min_flatness", "25"}}},
    {"TrialModel", {{}}},
    {"Log", {{"trials_per_write", "1"}, {"output_file", "tmp/me.txt"}}},
    {"Movie", {{"trials_per_write", "1"}, {"output_file", "tmp/me.xyz"}}},
    {"CheckEnergy", {{"trials_per_update", "1"}, {"tolerance", "1e-9"}}},
    {"Tune", {{}}},
    //{"Run", {{"num_trials", "1e2"}}},
  }}, true);
  auto mc2 = test_serialize_unique(*mc);
  mc2->run_num_trials(1e2);
}

}  // namespace feasst
