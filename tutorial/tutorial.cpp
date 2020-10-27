#include "feasst.h"

static feasst::ArgumentParse args("A canonical ensemble Metropolis Monte Carlo simulation of a bulk Lennard Jones fluid.", {
  {"--task", "SLURM job array index", "0"},
  {"--seed", "random number generator seed", "time"},
  {"--length", "cubic periodic boundary length", "8"},
  {"--num", "number of particles", "50"},
  {"--data", "LMP forcefield data file",
    feasst::install_dir() + "/forcefield/data.lj"},
  {"--beta", "inverse temperature", "1.2"},
  {"--trials", "number of Monte Carlo trials", feasst::str(1e6)}});

int main(int argc, char ** argv) {
  INFO(feasst::version());
  INFO("args: " << args.parse(argc, argv));
  const int num_trials = args.get_int("--trials");
  if (args.get_int("--task") > 0) {
    auto mc = feasst::MakeMonteCarlo("checkpoint.fst");
    INFO(mc->trials().num_attempts());
    mc->attempt(num_trials - mc->trials().num_attempts());
    return 0;
  }
  feasst::MonteCarlo mc;
  mc.set(feasst::MakeRandomMT19937({{"seed", args.get("--seed")}}));
  mc.add(feasst::Configuration(
    feasst::MakeDomain({{"cubic_box_length", args.get("--length")}}),
    {{"particle_type", args.get("--data")}}));
  mc.add(feasst::Potential(feasst::MakeLennardJones()));
  mc.add(feasst::Potential(feasst::MakeLongRangeCorrections()));
  mc.set(feasst::MakeThermoParams({{"beta", args.get("--beta")}}));
  mc.set(feasst::MakeMetropolis());
  mc.add(feasst::MakeTrialTranslate(
    {{"tunable_param", "2."}, {"tunable_target_acceptance", "0.2"}}));
  mc.add(feasst::MakeCheckEnergyAndTune(
   {{"steps_per", feasst::str(1e5)}, {"tolerance", "1e-8"}}));
  mc.set(feasst::MakeCheckpoint({{"file_name", "checkpoint.fst"},
                                 {"num_hours", "0.01"},
                                 {"num_hours_terminate", "0.03"}}));
  feasst::SeekNumParticles(args.get_int("--num"))
    .with_thermo_params({{"beta", "0.1"}, {"chemical_potential", "10"}})
    .with_metropolis()
    .with_trial_add()
    .run(&mc);
  mc.add(feasst::MakeLogAndMovie(
   {{"steps_per", feasst::str(1e5)}, {"file_name", "lj"}}));
  mc.attempt(num_trials);
}
