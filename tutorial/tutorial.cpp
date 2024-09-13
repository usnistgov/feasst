/**
 * This is an example usage of the C++ interface of FEASST.
 *
 * Usage:
 *   mkdir build; cd &_
 *   cmake ..
 *   make
 *   ./tutorial
 */

#include <iostream>
#include "feasst.h"

static feasst::ArgumentParse args(
  "Canonical ensemble Metropolis Monte Carlo simulation of Lennard Jones.", {
  {"--seed", "random number generator seed", "time"},
  {"--length", "cubic periodic boundary length", "8"},
  {"--num", "number of particles", "50"},
  {"--data", "LMP particle data file",
    feasst::install_dir() + "/particle/lj.fstprt"},
  {"--beta", "inverse temperature", "1.2"},
  {"--trials", "number of Monte Carlo trials", "1e6"}});

int main(int argc, char ** argv) {
  std::cout << "#FEASST version: " << feasst::version() << std::endl
            << args.parse(argc, argv) << std::endl;

  // high temperature gcmc to generate initial configuration
  auto mc = feasst::MakeMonteCarlo({{
    {"RandomMT19937", {{"seed", args.get("--seed")}}},
    {"Configuration", {{"cubic_side_length", args.get("--length")},
                       {"particle_type", args.get("--data")}}},
    {"Potential", {{"Model", "LennardJones"}}},
    {"Potential", {{"VisitModel", "LongRangeCorrections"}}},
    {"ThermoParams", {{"beta", "0.1"}, {"chemical_potential", "10"}}},
    {"Metropolis", {{}}},
    {"TrialTranslate", {{"tunable_param", "2."},
                        {"tunable_target_acceptance", "0.2"}}},
    {"TrialAdd", {{"particle_type", "0"}}},
    {"Run", {{"until_num_particles", "50"}}}
  }});

  ASSERT(mc->configuration().num_particles() == 50,
    "There should be 50 particles");
  ASSERT(mc->trials().num() == 2, "There should be two Trials");

  // nvt equilibration
  mc->begin({{
    {"RemoveTrial", {{"name", "TrialAdd"}}},
    {"ThermoParams", {{"beta", args.get("--beta")}}},
    {"CheckEnergy", {{"trials_per_update", "1e5"}, {"tolerance", "1e-8"}}},
    {"Tune", {{}}},
    {"Run", {{"num_trials", "1e5"}}}
  }});

  ASSERT(mc->trials().num() == 1, "There should be one Trial");

  // nvt production
  mc->begin({{
    {"Log", {{"trials_per_write", "1e5"}, {"output_file", "lj.csv"}}},
    {"Movie", {{"trials_per_write", "1e5"}, {"output_file", "lj.xyz"}}},
    {"Energy", {{"trials_per_write", "1e5"}, {"output_file", "lj_en.csv"}}},
    {"Run", {{"num_trials", args.get("--trials")}}}
  }});

  std::cout << "Ensemble average extensive potential energy: " << std::endl <<
    mc->analyze(2).accumulator().str() << std::endl;
}
