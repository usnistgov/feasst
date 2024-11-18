#include "utils/test/utils.h"
#include "monte_carlo/include/monte_carlo.h"
#include "chain/include/perturb_branch.h"

namespace feasst {

TEST(PerturbBranch, serialize) {
  PerturbBranch branch;
  PerturbBranch branch2 = test_serialize(branch);
}

TEST(PerturbBranch, tip4p) {
//  System sys;
//  sys.add(MakeConfiguration({{"cubic_side_length", "10"},
//    {"particle_type", "../plugin/charge/particle/tip4p.fstprt"},
//    {"add_particles_of_type0", "1"}}));
//  sys.add(MakePotential(MakeHardSphere(), MakeDontVisitModel()));
//  sys.precompute();
//  PerturbBranch perturb;
//  auto sel = MakeTrialSelectBranch({{"mobile_site", "1"}, {"mobile_site2", "2"}, {"anchor_site", "0"}, {"anchor_site2", "3"}});
//  sel.select_particle(0, sys.configuration());
}

TEST(MonteCarlo, tip4p) {
  auto mc = MakeMonteCarlo({{
    //{"RandomMT19937", {{"seed", "123"}}},
    {"Configuration", {{"particle_type0", "../plugin/charge/particle/tip4p.fstprt"},
      {"cubic_side_length", "20"}, {"add_particles_of_type0", "1"}}},
    {"Potential", {{"Model", "HardSphere"}, {"VisitModel", "DontVisitModel"}}},
    {"ThermoParams", {{"beta", "1"}, {"chemical_potential0", "-1"}}},
    {"Metropolis", {{}}},
    {"TrialGrowFile", {{"grow_file", "../plugin/chain/test/data/tip4p.grow"}}},
    {"TrialParticlePivot", {{"particle_type", "0"}}},
    {"TrialTranslate", {{}}},
//    {"Log", {{"trials_per_write", str(1e0)}, {"output_file", "tmp/tip4p.txt"}}},
//    {"Movie", {{"trials_per_write", str(1e0)}, {"output_file", "tmp/tip4p.xyz"}}},
//    {"CheckEnergy", {{"trials_per_update", str(1e0)}, {"tolerance", str(1e-9)}}},
//    {"Tune", {{}}},
    {"Run", {{"num_trials", "1e3"}}},
  }}, true);
}

}  // namespace feasst
