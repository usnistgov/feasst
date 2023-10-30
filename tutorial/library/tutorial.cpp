#include <iostream>
#include "feasst.h"

static feasst::ArgumentParse args("A canonical ensemble Metropolis Monte Carlo simulation of a bulk Lennard Jones fluid.", {
  {"--seed", "random number generator seed", "time"},
  {"--length", "cubic periodic boundary length", "8"},
  {"--num", "number of particles", "50"},
  {"--data", "LMP particle data file",
    feasst::install_dir() + "/particle/lj.fstprt"},
  {"--beta", "inverse temperature", "1.2"},
  {"--trials", "number of Monte Carlo trials", "1e6"}});

int main(int argc, char ** argv) {
  std::cout << feasst::version() << std::endl
            << args.parse(argc, argv) << std::endl;

  // high temperature gcmc to generate initial configuration
  feasst::MonteCarlo mc;
  mc.set(feasst::MakeRandomMT19937({{"seed", args.get("--seed")}}));
  mc.add(feasst::MakeConfiguration({{"cubic_side_length", args.get("--length")},
                                    {"particle_type", args.get("--data")}}));
  mc.add(feasst::MakePotential({{"Model", "LennardJones"}}));
  mc.add(feasst::MakePotential({{"VisitModel", "LongRangeCorrections"}}));
  mc.set(feasst::MakeThermoParams({{"beta", "0.1"}, {"chemical_potential", "10"}}));
  mc.set(feasst::MakeMetropolis());
  mc.add(feasst::MakeTrialTranslate(
    {{"tunable_param", "2."}, {"tunable_target_acceptance", "0.2"}}));
  mc.add(feasst::MakeTrialAdd({{"particle_type", "0"}}));
  mc.run(feasst::MakeRun({{"until_num_particles", "50"}}));

  // nvt equilibration
  mc.run(feasst::MakeRemoveTrial({{"name", "TrialAdd"}}));
  mc.set(feasst::MakeThermoParams({{"beta", args.get("--beta")}}));
  mc.add(feasst::MakeCheckEnergy({{"trials_per_update", "1e5"}, {"tolerance", "1e-8"}}));
  mc.add(feasst::MakeTune());
  mc.run(feasst::MakeRun({{"num_trials", "1e5"}}));

  // nvt production
  mc.add(feasst::MakeLog({{"trials_per_write", "1e5"}, {"output_file", "lj.csv"}}));
  mc.add(feasst::MakeMovie({{"trials_per_write", "1e5"}, {"output_file", "lj.xyz"}}));
  mc.run(feasst::MakeRun({{"num_trials", args.get("--trials")}}));
}
