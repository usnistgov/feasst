#include <iostream>
#include "feasst.h"

static feasst::ArgumentParse args("A canonical ensemble Metropolis Monte Carlo simulation of a bulk Lennard Jones fluid.", {
  {"--task", "SLURM job array index", "0"},
  {"--seed", "random number generator seed", "time"},
  {"--density", "number density, N/V", "0.776"},
  {"--num", "number of particles", "500"},
  {"--data", "FEASST particle file",
    feasst::install_dir() + "/particle/lj.fstprt"},
  {"--beta", "inverse temperature", feasst::str(1./0.9)},
  {"--equilibration_trials", "number of equilibration trials", feasst::str(5e7)},
  {"--trials", "total number of trials (equilibration and production)", feasst::str(3e8)},
  {"--num_hours", "number of hours before restart", "1"}});

int main(int argc, char ** argv) {
  std::cout << feasst::version() << std::endl
            << args.parse(argc, argv) << std::endl;
  const int num_trials = args.get_int("--trials");
  if (args.get_int("--task") > 0) {
    auto mc = feasst::MakeMonteCarlo("checkpoint.fst");
    std::cout << mc->trials().num_attempts() << std::endl;
    mc->attempt(num_trials - mc->trials().num_attempts());
    return 0;
  }
  feasst::MonteCarlo mc;
  mc.set(feasst::MakeRandomMT19937({{"seed", args.get("--seed")}}));
  mc.add(feasst::Configuration(
    feasst::MakeDomain({{"cubic_side_length", feasst::str(std::pow(args.get_double("--num")/args.get_double("--density"), 1./3.))}}),
    {{"particle_type", args.get("--data")}}));
  mc.add(feasst::MakePotential(feasst::MakeLennardJones()));
  mc.add(feasst::MakePotential(feasst::MakeLongRangeCorrections()));
  mc.set(feasst::MakeThermoParams({{"beta", args.get("--beta")}}));
  mc.set(feasst::MakeMetropolis());
  mc.add(feasst::MakeTrialTranslate(
    {{"tunable_param", "0.2"}, {"tunable_target_acceptance", "0.2"}}));
  mc.add(feasst::MakeTrialAdd({{"particle_type", "0"}}));
  mc.run(feasst::MakeRun({{"until_num_particles", args.get("--num")}}));
  mc.run(feasst::MakeRemoveTrial({{"name", "TrialAdd"}}));
  mc.add(feasst::MakeCheckEnergy({{"trials_per", feasst::str(1e5)}, {"tolerance", "1e-8"}}));
  mc.add(feasst::MakeTune());
  mc.set(feasst::MakeCheckpoint({{"file_name", "checkpoint.fst"},
                                 {"num_hours", feasst::str(0.95*args.get_double("--num_hours"))},
                                 {"num_hours_terminate", feasst::str(0.95*args.get_double("--num_hours"))}}));
  mc.add(feasst::MakeLogAndMovie(
    {{"trials_per", feasst::str(1e5)}, {"file_name", "lj"}}));
  mc.add(feasst::MakeIncrementPhase({{"num_trials", args.get("--equilibration_trials")}}));
  mc.add(feasst::MakeEnergy(
    {{"trials_per_write", feasst::str(1e5)}, {"file_name", "en_lj.txt"}, {"start_after_phase", "0"}}));
  mc.attempt(num_trials);
}
