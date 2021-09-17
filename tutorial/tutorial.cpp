#include <iostream>
#include "feasst.h"

static feasst::ArgumentParse args("A canonical ensemble Metropolis Monte Carlo simulation of a bulk Lennard Jones fluid.", {
  {"--seed", "random number generator seed", "time"},
  {"--length", "cubic periodic boundary length", "8"},
  {"--num", "number of particles", "50"},
  {"--data", "LMP forcefield data file",
    feasst::install_dir() + "/forcefield/lj.fstprt"},
  {"--beta", "inverse temperature", "1.2"},
  {"--trials", "number of Monte Carlo trials", "1e6"}});

int main(int argc, char ** argv) {
  std::cout << feasst::version() << std::endl
            << args.parse(argc, argv) << std::endl;

  // high temperature gcmc to generate initial configuration
  feasst::MonteCarlo mc;
  mc.set(feasst::MakeRandomMT19937({{"seed", args.get("--seed")}}));
  mc.add(feasst::MakeConfiguration({{"cubic_box_length", args.get("--length")},
                                    {"particle_type", args.get("--data")}}));
  mc.add(feasst::MakePotential({{"Model", "LennardJones"}}));
  mc.add(feasst::MakePotential({{"VisitModel", "LongRangeCorrections"}}));
  mc.set(feasst::MakeThermoParams({{"beta", "0.1"}, {"chemical_potential", "10"}}));
  mc.set(feasst::MakeMetropolis());
  mc.add(feasst::MakeTrialTranslate(
    {{"tunable_param", "2."}, {"tunable_target_acceptance", "0.2"}}));
  mc.add(feasst::MakeTrialAdd({{"particle_type", "0"}}));
  mc.perform(feasst::MakeRun({{"until_num_particles", "50"}}));

  // nvt equilibration
  mc.perform(feasst::MakeRemoveTrial({{"name", "TrialAdd"}}));
  mc.set(feasst::MakeThermoParams({{"beta", args.get("--beta")}}));
  mc.add(feasst::MakeCheckEnergyAndTune(
   {{"steps_per", "1e5"}, {"tolerance", "1e-8"}}));
  mc.perform(feasst::MakeRun({{"num_attempts", "1e5"}}));

  // nvt production
  mc.add(feasst::MakeLogAndMovie(
   {{"steps_per", "1e5"}, {"file_name", "lj"}}));
  mc.perform(feasst::MakeRun({{"num_attempts", args.get("--trials")}}));
}
