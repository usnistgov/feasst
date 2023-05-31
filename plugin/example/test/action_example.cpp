#include <cmath>
#include "utils/test/utils.h"
#include "monte_carlo/include/monte_carlo.h"
#include "example/include/action_example.h"

namespace feasst {

TEST(ActionExample, serialize) {
  auto action = MakeActionExample();
  auto action2 = test_serialize<ActionExample, Action>(*action);
}

TEST(ActionExample, fh) {
  auto mc = MakeMonteCarlo({{
    {"Configuration", {{"particle_type0", "../forcefield/atom.fstprt"},
                       {"cubic_box_length", "8"}}},
    {"Potential", {{"Model", "IdealGas"}}},
    {"ThermoParams", {{"beta", "1"}, {"chemical_potential0", "1"}}},
    {"FlatHistogram", {{"Macrostate", "MacrostateNumParticles"}, {"width", "1"}, {"max", "4"}, {"min", "0"}, {"Bias", "TransitionMatrix"}, {"min_sweeps", "2"}}},
    {"TrialTransfer", {{"particle_type", "0"}}},
    {"CriteriaUpdater", {{"trials_per_update", "1e2"}}},
    {"CriteriaWriter", {{"trials_per_write", "1e2"}, {"file_name", "tmp/action_ex_crit.csv"}}},
    {"Energy", {{"file_name", "tmp/en.csv"}, {"trials_per_write", "1e2"}, {"multistate", "true"}}},
    {"PairDistribution", {{"file_name", "tmp/grig.csv"}, {"trials_per_write", "1e2"}}},
    {"Run", {{"until_criteria_complete", "true"}}},
    {"ActionExample", {{"analyze_name", "CriteriaWriter"}}},
    {"ActionExample", {{"analyze_name", "Energy"}}},
    {"ActionExample", {{"modify_name", "PairDistribution"}}},
  }});
}

}  // namespace feasst
