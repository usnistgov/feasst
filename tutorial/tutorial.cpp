#include "feasst.h"

feasst::ArgumentParse args(
"A canonical ensemble Metropolis Monte Carlo simulation of a bulk Lennard "
"Jones fluid.\n\n"

"options:\n"
"--help   : print this documentation to screen.\n"
"--task   : SLURM job array index (default: 0).\n"
"--seed   : random number generator seed (default: time).\n"
"--length : cubic periodic boundary length (default: 8).\n"
"--num    : number of particles (default: 50).\n"
"--data   : LMP forcefield data file (default: feasst/forcefield/data.lj).\n"
"--beta   : inverse temperature (default: 1.2).\n"
"--trials : number of Monte Carlo trials (default: 1e6).\n"
);

int main(int argc, char ** argv) {
  INFO(feasst::version());
  INFO("args: " << args.parse(argc, argv));
  const int num_trials = args.get_int("--trials", 1e6);
  if (args.get_int("--task", 0) > 0) {
    auto mc = feasst::MakeMonteCarlo("checkpoint.fst");
    INFO(mc->trials().num_attempts());
    mc->attempt(num_trials - mc->trials().num_attempts());
    return 0;
  }
  feasst::MonteCarlo mc;
  mc.set(feasst::MakeRandomMT19937({{"seed", args.get("--seed", "time")}}));
  mc.add(feasst::Configuration(
    feasst::MakeDomain({{"cubic_box_length", args.get("--length", "8")}}),
    {{"particle_type",
      args.get("--data", feasst::install_dir() + "/forcefield/data.lj")}}));
  mc.add(feasst::Potential(feasst::MakeLennardJones()));
  mc.add(feasst::Potential(feasst::MakeLongRangeCorrections()));
  mc.add(feasst::MakeMetropolis({{"beta", args.get("--beta", "1.2")}}));
  mc.add(feasst::MakeTrialTranslate(
    {{"tunable_param", "2."}, {"tunable_target_acceptance", "0.2"}}));
  mc.add(feasst::MakeCheckEnergyAndTune(
   {{"steps_per", feasst::str(1e5)}, {"tolerance", "1e-8"}}));
  mc.set(feasst::MakeCheckpoint({{"file_name", "checkpoint.fst"},
                                 {"num_hours", "0.01"},
                                 {"num_hours_terminate", "0.03"}}));
  feasst::SeekNumParticles(args.get_int("--num", 50))
    .with_metropolis({{"beta", "0.1"}, {"chemical_potential", "10"}})
    .with_trial_add()
    .run(&mc);
  mc.add(feasst::MakeLogAndMovie(
   {{"steps_per", feasst::str(1e5)}, {"file_name", "lj"}}));
  mc.attempt(num_trials);
}
