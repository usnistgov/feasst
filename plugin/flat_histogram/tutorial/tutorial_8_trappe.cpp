#include <fstream>
#include "feasst.h"

feasst::ArgumentParse args(
"A grand canonical ensemble flat histogram Monte Carlo simulation of a bulk fluid.\n\n"

"options:\n"
"--help/-h          : print this documentation to screen.\n"
"--task             : SLURM job array index (default=0).\n"
"--num_procs        : number of processors (default=12).\n"
"--num_hours        : number of hours before restart (default=5*24).\n"
"--dccb_cutoff      : cutoff for dual cut configurational bias (default=4).\n"
"--dccb_begin       : number of molecules before using DCCB (default=300).\n"
"--lx               : box length in x (default=33.0).\n"
"--ly               : box length in y (default=33.0).\n"
"--lz               : box length in z (default=33.0).\n"
"--max_particles    : maximum number of particles (default=400).\n"
"--temperature      : temperature in Kelvin (default=263.15).\n"
"--particle0        : particle type 0.\n"
"--particle1        : particle type 1, not used if empty (default).\n"
"--collect_flatness : number of WL flatness to begin collection (default=18).\n"
"--min_flatness     : number of WL flatness to switch to TM (default=22).\n"
"--beta_mu          : baseline chemical potential of each species (default=-7).\n"
"--delta_betamu_0   : delta_betamu_0 (default=0.)\n"
"--cyl_radius       : radius of cylindrical confinement, not used if < 0 (default: -1).\n"
"--cyl_cutoff       : square well interaction distance from cylinder to point (default: 6).\n"
"--cyl_epsilon      : square well interaction strength (default: 500).\n"
);

std::shared_ptr<feasst::MonteCarlo> mc(const int thread, const int mn, const int mx) {
  const std::string steps_per = feasst::str(1e4);
  auto mc = feasst::MakeMonteCarlo();
  feasst::argtype domain_args = {
    {"side_length0", args.get("--lx", "33")},
    {"side_length1", args.get("--ly", "33")},
    {"side_length2", args.get("--lz", "33")}};
  const int dccb_begin = args.get_int("--dccb_begin", 300);
  std::string ref("-1");
  std::string num("1");
  const int dccb_cutoff = args.get_int("--dccb_cutoff", 4);
  if (mx > dccb_begin) {
    domain_args.insert({"init_cells", feasst::str(dccb_cutoff)});
    ref = "0";
    num = "4";
  }
  const double beta = 1./args.get_double("--temperature", 263.15);
  const double beta_mu = args.get_double("--beta_mu", -7.);
  feasst::argtype config_args = {{"particle_type0", args.get("--particle0")}};
  feasst::argtype criteria_args = {{"beta", feasst::str(beta)},
    {"chemical_potential0", feasst::str((beta_mu+args.get_double("--delta_betamu_0", 0))/beta)}};
  const std::string part1 = args.get("--particle1");
  if (!part1.empty()) {
    config_args.insert({"particle_type1", part1});
    criteria_args.insert({"chemical_potential1", feasst::str(beta_mu/beta)});
  }
  mc->add(feasst::Configuration(feasst::MakeDomain(domain_args), config_args));
  mc->add(feasst::Potential(feasst::MakeLennardJones()));
  mc->add(feasst::Potential(feasst::MakeLongRangeCorrections()));
  const double cyl_radius = args.get_double("--cyl_radius", -1.);
  if (cyl_radius > 0) {
    feasst::Potential cylinder(feasst::MakeModelSquareWellShape(feasst::MakeCylinder(
      {{"radius", feasst::str(cyl_radius)}},
      feasst::Position({{"x", "0"}, {"y", "0"}, {"z", "0"}}),
      feasst::Position({{"x", "0"}, {"y", "0"}, {"z", "1"}}))));
    cylinder.set_model_params(mc->configuration());
    for (int site_type = 0; site_type < mc->configuration().num_site_types(); ++site_type) {
      cylinder.set_model_param("cutoff", site_type, args.get_double("--cyl_cutoff", 10));
      cylinder.set_model_param("epsilon", site_type, args.get_double("--cyl_epsilon", 500));
    }
    mc->add(cylinder);
  }
  if (mx > dccb_begin) {
    feasst::Potential reference(feasst::MakeLennardJones());
    if (mc->configuration().domain().num_cells() > 0) {
      reference = feasst::Potential(feasst::MakeLennardJones(), feasst::MakeVisitModelCell());
    }
    reference.set_model_params(mc->configuration());
    for (int site_type = 0; site_type < mc->configuration().num_site_types(); ++site_type) {
      reference.set_model_param("cutoff", site_type, dccb_cutoff);
    }
    mc->add_to_reference(reference);
  }
  mc->add(feasst::MakeFlatHistogram(
    feasst::MakeMacrostateNumParticles(
      feasst::Histogram({{"width", "1"}, {"max", feasst::str(mx)}, {"min", feasst::str(mn)}})),
    // feasst::MakeTransitionMatrix(feasst::args({"min_sweeps": feasst::str(args.sweeps)})),
    feasst::MakeWLTM({
      {"collect_flatness", args.get("--collect_flatness", "18")},
      {"min_flatness", args.get("--min_flatness", "22")},
      {"min_sweeps", "1000"}}),
    criteria_args));
  for (int particle_type = 0; particle_type < mc->configuration().num_particle_types(); ++particle_type) {
    mc->add(feasst::MakeTrialTranslate({
      {"particle_type", feasst::str(particle_type)},
      {"weight", "1."},
      {"tunable_param", "1."},
      {"reference_index", ref},
      {"num_steps", num}}));
    mc->add(feasst::MakeTrialRotate({
      {"particle_type", feasst::str(particle_type)},
      {"weight", "1."},
      {"tunable_param", "1."},
      {"reference_index", ref},
      {"num_steps", num}}));
    mc->add(feasst::MakeTrialTransfer({
      {"particle_type", feasst::str(particle_type)},
      {"weight", "4"},
      {"reference_index", ref},
      {"num_steps", num}}));
  }
  mc->add(feasst::MakeCheckEnergy({{"steps_per", steps_per}, {"tolerance", "0.0001"}}));
  mc->add(feasst::MakeTuner({{"steps_per", steps_per}, {"stop_after_phase", "0"}}));
  // mc->add(feasst::MakeLogAndMovie({{"steps_per", steps_per},
  //                                  {"file_name", "clones" + feasst::str(thread)},
  //                                  {"file_name_append_phase", "True"}}));
  mc->add(feasst::MakeEnergy({
    {"file_name", "en" + feasst::str(thread) + ".txt"},
    {"file_name_append_phase", "true"},
    {"start_after_phase", "0"},
    {"steps_per_write", steps_per},
    {"steps_per_update", "1"},
    {"multistate", "True"}}));
  if (mc->configuration().num_particle_types() > 1) {
    mc->add(feasst::MakeNumParticles({
      {"particle_type", "0"},
      {"file_name", "num" + feasst::str(thread) + ".txt"},
      {"file_name_append_phase", "true"},
      {"start_after_phase", "0"},
      {"steps_per_write", steps_per},
      {"steps_per_update", "1"},
      {"multistate", "True"}}));
  }
  mc->add(feasst::MakeExtensiveMoments({
    {"steps_per_write", steps_per},
    {"file_name", "extmom"+feasst::str(thread) + ".txt"},
    {"file_name_append_phase", "true"},
    {"start_after_phase", "0"}, // do not update until equilibration complete (phase 1)
    {"max_order", "3"},
    {"multistate", "True"}}));
  mc->add(feasst::MakeCriteriaUpdater({{"steps_per", steps_per}}));
  mc->add(feasst::MakeCriteriaWriter({{"steps_per", steps_per},
                                     {"file_name", "clones" + feasst::str(thread) + "_crit.txt"},
                                     {"file_name_append_phase", "True"}}));
  mc->set(feasst::MakeCheckpoint({{"file_name", "checkpoint" + feasst::str(thread) + ".fst"},
                                 {"num_hours_terminate", feasst::str(0.9*args.get_int("--num_procs", 12)*args.get_double("--num_hours", 5*24))}}));
  mc->add(feasst::MakeAnalyzeRigidBonds({{"steps_per", steps_per}}));
  return mc;
}

int main(int argc, char ** argv) {
  INFO(feasst::version());
  INFO("args: " << args.parse(argc, argv));
  const int num_procs = args.get_int("--num_procs", 12);
  auto windows = feasst::WindowExponential({
    {"alpha", "1.75"},
    {"num", feasst::str(num_procs)},
    {"maximum", args.get("--max_particles", "400")},
    {"extra_overlap", "2"}}).boundaries();
  INFO(feasst::feasst_str(windows));
  auto clones = feasst::MakeClones();
  if (args.get_int("--task", 0) == 0) {
    for (int proc = 0; proc < static_cast<int>(windows.size()); ++proc) {
      clones->add(mc(proc, windows[proc][0], windows[proc][1]));
    }
    clones->set(feasst::MakeCheckpoint({{"file_name", "checkpoint.fst"}}));
  } else {
    clones = feasst::MakeClones("checkpoint", num_procs);
  }
  clones->initialize_and_run_until_complete({{"ln_prob_file", "ln_prob.txt"}});
  INFO(feasst::feasst_str(clones->ln_prob().values()));
  std::ofstream file("clones.fst");
  file << clones->serialize();
  file.close();
}
