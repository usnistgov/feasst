#include "utils/test/utils.h"
#include "math/include/accumulator.h"
#include "monte_carlo/include/analyze.h"
#include "monte_carlo/include/monte_carlo.h"
#include "flat_histogram/include/flat_histogram.h"
#include "flat_histogram/include/ln_probability.h"
#include "flat_histogram/include/bias.h"
#include "chain/include/trial_reptate_unopt.h"

namespace feasst {

TEST(TrialReptateUnopt, serialize) {
  auto trial = std::make_shared<TrialReptateUnoptHalf>();
  Trial trial2 = test_serialize(*trial);
}

//TEST(MonteCarlo, reptate) {
//  auto mc = MakeMonteCarlo({{
//    {"Configuration", {{"cubic_side_length", "30"}, {"particle_type1", "../particle/n-octane.fstprt"},
//                       {"add_particles_of_type1", "1"}, {"particle_type0", "../particle/atom.fstprt"},
//                       {"sigma0", "4"}}},
//    {"Potential", {{"Model", "LennardJones"}}},
//    {"Potential", {{"Model", "LennardJones"}, {"VisitModel", "VisitModelIntra"}, {"intra_cut", "3"}}},
//    {"ThermoParams", {{"beta", "1"}, {"chemical_potential0", "1."}, {"chemical_potential1", "1."}}},
//    {"Metropolis", {{}}},
//    {"TrialReptateUnopt", {{"particle_type", "1"}}},
//    {"TrialTranslate", {{"particle_type", "0"}, {"tunable_param", "1"}}},
//    {"TrialAdd", {{"particle_type", "0"}}},
//    {"Run", {{"until_num_particles", "40"}}},
//    {"Remove", {{"name", "TrialAdd"}}},
//    {"Log", {{"trials_per_write", str(1e0)}, {"output_file", "tmp/rept.txt"}}},
//    {"Movie", {{"trials_per_write", str(1e0)}, {"output_file", "tmp/rept.xyz"}}},
//    {"CheckEnergy", {{"trials_per_update", str(1e0)}, {"tolerance", str(1e-9)}}},
//    {"Tune", {{}}},
//    {"Run", {{"num_trials", "1e2"}}},
//  }}, true);
//}
//
//TEST(MonteCarlo, reptate_fh) {
//  const std::string tpc = "1e4";
//  auto mc = MakeMonteCarlo({{
//    {"RandomMT19937", {{"seed", "123"}}},
//    {"Configuration", {{"cubic_side_length", "8"}, {"particle_type0", "../plugin/chain/test/data/trimer_harmonic.fstprt"}, {"add_particles_of_type0", "1"}}},
//    //{"Configuration", {{"cubic_side_length", "45"}, {"particle_type0", "../particle/n-octane.fstprt"}, {"add_particles_of_type0", "1"}}},
//    {"Potential", {{"Model", "LennardJones"}}},
//    //{"Potential", {{"Model", "LennardJones"}, {"VisitModel", "VisitModelIntraMap"}, {"exclude_bonds", "true"}, {"exclude_angles", "true"}, {"exclude_dihedrals", "true"}}},
//    {"RefPotential", {{"VisitModel", "DontVisitModel"}}},
//    {"ThermoParams", {{"beta", "1"}}},
//    //{"ThermoParams", {{"beta", "1"}, {"chemical_potential0", "0"}}},
//    //{"ThermoParams", {{"beta", "0.3436353001220744"}, {"chemical_potential0", "-17.460371498121805"}}},
//    {"Metropolis", {{}}},
//    {"TrialReptateUnopt", {{"particle_type", "0"}, {"reference_index", "0"}}},
//    {"TrialTranslate", {{"weight", "0.01"}}}, // to wrap particle in box
//    //{"TrialGrowFile", {{"grow_file", "../plugin/chain/test/data/trappe_grow_canonical.txt"}}},
//    //{"TrialGrowFile", {{"grow_file", "../plugin/chain/test/data/trimer_grow_canonical.txt"}}},
//    {"Energy", {{"trials_per_write", tpc}, {"output_file", "tmp/rept_en.txt"}}},
//    {"Log", {{"trials_per_write", tpc}, {"output_file", "tmp/rept.txt"}}},
//    {"Movie", {{"trials_per_write", tpc}, {"output_file", "tmp/rept.xyz"}}},
//    {"CheckEnergy", {{"trials_per_update", tpc}, {"tolerance", str(1e-9)}}},
//    {"Tune", {{}}},
//    {"Run", {{"num_trials", "1e1"}}},
//    //{"Run", {{"num_trials", "1e6"}}},
//    //{"Run", {{"until", "complete"}}},
//  }}, true);
//  const Analyze * en = mc->analyzers()[0].get();
//  EXPECT_NEAR(en->accumulator().average(), 0.57168889603803608, 5*0.0014774057248454694);
//  INFO(mc->trial(0).acceptance() << " " << mc->trial(0).num_success());
//  INFO("add reptation to equipartition tests, although, on average, bond lengths are slightly longer than equilibrium length");
//}

}  // namespace feasst
