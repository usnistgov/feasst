#include <cmath>
#include "utils/test/utils.h"
#include "monte_carlo/include/monte_carlo.h"
#include "steppers/include/write_stepper.h"

namespace feasst {

TEST(WriteStepper, serialize) {
  auto action = MakeWriteStepper();
  auto action2 = test_serialize<WriteStepper, Action>(*action);
}

TEST(WriteStepper, fh) {
  auto mc = MakeMonteCarlo({{
    {"Configuration", {{"particle_type0", "../particle/atom.fstprt"},
                       {"cubic_side_length", "8"}}},
    {"Potential", {{"Model", "IdealGas"}}},
    {"ThermoParams", {{"beta", "1"}, {"chemical_potential0", "1"}}},
    {"FlatHistogram", {{"Macrostate", "MacrostateNumParticles"}, {"width", "1"}, {"max", "4"}, {"min", "0"}, {"Bias", "TransitionMatrix"}, {"min_sweeps", "2"}}},
    {"TrialTransfer", {{"particle_type", "0"}}},
    {"CriteriaUpdater", {{"trials_per_update", "1e2"}}},
    {"CriteriaWriter", {{"trials_per_write", "1e2"}, {"output_file", "tmp/action_ex_crit.csv"}}},
    {"Energy", {{"output_file", "tmp/en.csv"}, {"trials_per_write", "1e2"}, {"multistate", "true"}}},
    {"PairDistribution", {{"output_file", "tmp/grig.csv"}, {"trials_per_write", "1e2"}}},
    {"Run", {{"until_criteria_complete", "true"}}},
    {"WriteStepper", {{"analyze_name", "CriteriaWriter"}}},
    {"WriteStepper", {{"analyze_name", "Energy"}}},
    {"WriteStepper", {{"modify_name", "PairDistribution"}}},
  }}, true);
}

}  // namespace feasst
