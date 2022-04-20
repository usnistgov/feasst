#include <unistd.h>
#include <iostream>
#include <sstream>
#include "feasst.h"

using namespace feasst;

System lj_system(const double cutoff) {
  System system;
  std::stringstream particle;
  particle << install_dir() << "/forcefield/lj.fstprt";
  Configuration config(
      MakeDomain({{"cubic_box_length", str(2.*cutoff)}}),
      {{"particle_type", particle.str()}});
  config.set_model_param("cutoff", 0, cutoff);
  system.add(config);
  system.add(MakePotential(MakeLennardJones()));
  return system;
}

int main(int argc, char** argv) {
  int pipe = true;
  double cutoff = 3.;
  double target_prob = -1;
  double max_move = -1;

  // parse command-line arguments using getopt
  { int index, c; opterr = 0;
    while ((c = getopt(argc, argv, "p:c:t:m:")) != -1) {
      switch (c) {
        case 'p': pipe = atoi(optarg); break;
        case 'c': cutoff = atof(optarg); break;
        case 't': target_prob = atof(optarg); break;
        case 'm': max_move = atof(optarg); break;
	      case '?':
          if (isprint(optopt))
            fprintf(stderr, "Unknown option `-%c'.\n", optopt);
          else
            fprintf(stderr, "Unknown option character `\\x%x'.\n", optopt);
          return 1;
        default: abort();
      }
    }
    std::cout << "# -p pipe " << pipe << " -c cutoff " << cutoff << " -t target_prob " << target_prob << " -m max_move " << max_move << std::endl;
    for (index = optind; index < argc; index++) printf("Non-option argument %s\n", argv[index]);
  }  // GETOPT

  const double density=0.85;
  const double temperature=0.88;
  const int trials_per = 1e6;
  const double box_length = 2*cutoff;
  std::stringstream file_app;
  file_app << "_a" << target_prob << "_r" << cutoff;
  const int num_particles = density*pow(box_length, 3);

  auto monte_carlo = MakePrefetch({{"trials_per_check", "10000000"}});
  monte_carlo->activate_prefetch(false);
  monte_carlo->set(lj_system(cutoff));
  monte_carlo->set(MakeMetropolis({
    {"beta", str(1./temperature)},
    {"chemical_potential", "1."}}));
  monte_carlo->add(MakeTrialTranslate({
    {"weight", "1."},
    {"tunable_param", "2."},
    {"tunable_param", str(max_move)},
    {"tunable_target_acceptance", str(target_prob)},
    {"tunable_percent_change", "0.01"}}));
  monte_carlo->seek_num_particles(num_particles);
  std::stringstream logfile;
  logfile << "log" << file_app.str() << ".txt";
  monte_carlo->add(MakeLog({
    {"trials_per", str(trials_per)},
    {"file_name", logfile.str()},
    {"clear_file", "true"}}));
  monte_carlo->add(MakeCheckEnergy({
    {"trials_per", str(trials_per)},
    {"tolerance", str(1e-8)}}));
  monte_carlo->add(MakeTune({
    {"trials_per", str(trials_per)}}));

  // equilibrate
  monte_carlo->attempt(int(1e7));

  if (pipe == 1) monte_carlo->activate_prefetch(true);
  std::stringstream msdfile;
  msdfile << "msd" << file_app.str() << ".txt";
  monte_carlo->add(MakeMeanSquaredDisplacement({
    {"trials_per_update", "10000"},
    {"updates_per_origin", "1000"},
    {"file_name", msdfile.str()},
    {"trials_per_write", str(int(1e5))}}));

  std::stringstream cpufile;
  cpufile << "cpu" << file_app.str() << ".txt";
  monte_carlo->add(MakeCPUTime({
    {"trials_per_update", str(trials_per)},
    {"trials_per_write", str(trials_per)},
    {"file_name", cpufile.str()}}));

  monte_carlo->attempt(int(1e8));
  return 0;
}
