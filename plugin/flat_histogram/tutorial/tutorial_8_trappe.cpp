#include <assert.h>
#include <fstream>
#include "feasst.h"

feasst::ArgumentParse args("A grand canonical ensemble flat histogram Monte Carlo simulation of a bulk fluid", {
  {"--task", "SLURM job array index", "0"},
  {"--num_procs", "number of processors", "12"},
  {"--num_hours", "number of hours before restart", feasst::str(5*24)},
  {"--dccb_cutoff", "cutoff for dual cut configurational bias", "4"},
  {"--dccb_begin", "number of molecules before using DCCB", "300"},
  {"--lx", "box length in x", "33.0"},
  {"--ly", "box length in y", "33.0"},
  {"--lz", "box length in z", "33.0"},
  {"--min_particles", "minimum number of particles", "0"},
  {"--max_particles", "maximum number of particles", "400"},
  {"--temperature", "temperature in Kelvin", "263.15"},
  {"--particle0", "particle type 0"},
  {"--particle1", "particle type 1, not used if empty (the default)"},
  {"--collect_flatness", "number of WL flatness to begin collection", "18"},
  {"--min_flatness", "number of WL flatness to switch to TM", "22"},
  {"--beta_mu", "baseline chemical potential of each species", "-7"},
  {"--delta_beta_mu1", "beta_mu1 = beta_mu0 + delta_beta_mu1", "0."},
  {"--cyl_radius", "radius of cylindrical confinement, not used if < 0", "-1"},
  {"--cyl_cutoff", "square well interaction distance from cylinder to point", "10"},
  {"--cyl_epsilon", "square well interaction strength", "500"},
  {"--file_xyz", "file name of xyz positions of initial configuration. "
    "Assume all particle0, and do not use if empty (the default)"}});

std::shared_ptr<feasst::MonteCarlo> mc(const int thread, const int mn, const int mx) {
  const std::string steps_per = feasst::str(1e5);
  auto mc = feasst::MakeMonteCarlo();
  feasst::argtype domain_args = {
    {"side_length0", args.get("--lx")},
    {"side_length1", args.get("--ly")},
    {"side_length2", args.get("--lz")}};
  const int dccb_begin = args.get_int("--dccb_begin");
  std::string ref("-1");
  std::string num("1");
  const int dccb_cutoff = args.get_int("--dccb_cutoff");
  if (mx > dccb_begin) {
    domain_args.insert({"init_cells", feasst::str(dccb_cutoff)});
    ref = "0";
    num = "4";
  }
  const double beta = 1./args.get_double("--temperature");
  const double beta_mu = args.get_double("--beta_mu");
  feasst::argtype config_args = {{"particle_type0", args.get("--particle0")}};
  feasst::argtype criteria_args = {{"beta", feasst::str(beta)},
    {"chemical_potential0", feasst::str(beta_mu/beta)}};
  if (args.option_given("--particle1")) {
    config_args.insert({"particle_type1", args.get("--particle1")});
    criteria_args.insert({"chemical_potential1", feasst::str((beta_mu + args.get_double("--delta_beta_mu1"))/beta)});
  }
  if (args.option_given("--file_xyz")) {
    assert(thread == 0); // cannot initialize multiple threads to same configuration
    feasst::Configuration config(feasst::MakeDomain(domain_args), config_args);
    feasst::FileXYZ().load(args.get("--file_xyz"), &config);
    mc->add(config);
  } else {
    mc->add(feasst::Configuration(feasst::MakeDomain(domain_args), config_args));
  }
  mc->add(feasst::Potential(feasst::MakeLennardJones()));
  mc->add(feasst::Potential(feasst::MakeLongRangeCorrections()));
  const double cyl_radius = args.get_double("--cyl_radius");
  if (cyl_radius > 0) {
    feasst::Potential cylinder(feasst::MakeModelSquareWellShape(feasst::MakeCylinder(
      {{"radius", feasst::str(cyl_radius)}},
      feasst::Position({{"x", "0"}, {"y", "0"}, {"z", "0"}}),
      feasst::Position({{"x", "0"}, {"y", "0"}, {"z", "1"}}))));
    cylinder.set_model_params(mc->configuration());
    for (int site_type = 0; site_type < mc->configuration().num_site_types(); ++site_type) {
      cylinder.set_model_param("cutoff", site_type, args.get_double("--cyl_cutoff"));
      cylinder.set_model_param("epsilon", site_type, args.get_double("--cyl_epsilon"));
    }
    mc->add(cylinder);
  }
  if (mx > dccb_begin) {
    feasst::Potential reference(feasst::MakeLennardJones());
    if (mc->configuration().domain().num_cells() > 0) {
      reference = feasst::Potential(feasst::MakeLennardJones(),
                                    feasst::MakeVisitModelCell({{"min_length", feasst::str(dccb_cutoff)}}));
    }
    reference.set_model_params(mc->configuration());
    for (int site_type = 0; site_type < mc->configuration().num_site_types(); ++site_type) {
      reference.set_model_param("cutoff", site_type, dccb_cutoff);
    }
    mc->add_to_reference(reference);
  }
  mc->set(feasst::MakeThermoParams(criteria_args));
  mc->set(feasst::MakeFlatHistogram(
    feasst::MakeMacrostateNumParticles(
      feasst::Histogram({{"width", "1"}, {"max", feasst::str(mx)}, {"min", feasst::str(mn)}})),
    // feasst::MakeTransitionMatrix(feasst::args({"min_sweeps": feasst::str(args.sweeps)})),
    feasst::MakeWLTM({
      {"collect_flatness", args.get("--collect_flatness")},
      {"min_flatness", args.get("--min_flatness")},
      {"min_sweeps", "1000"}})));
  std::cout << "Initial energy: " << MAX_PRECISION << mc->criteria().current_energy() << std::endl;
  for (int particle_type = 0; particle_type < mc->configuration().num_particle_types(); ++particle_type) {
    mc->add(feasst::MakeTrialTranslate({{"particle_type", feasst::str(particle_type)}, {"weight", "1."},
      {"tunable_param", "1."}, {"reference_index", ref}, {"num_steps", num}}));
    if (mx > dccb_begin && mc->configuration().particle_type(particle_type).num_sites() == 2) {
      mc->add(feasst::MakeTrialGrow({
        {{"transfer", "true"},  // {"regrow", "true"},  # regrow isn't very efficient
         {"particle_type", feasst::str(particle_type)}, {"site", "0"}, {"weight", "4"}},
        {{"bond", "true"}, {"mobile_site", "1"}, {"anchor_site", "0"}}},
        {{"reference_index", ref}, {"num_steps", num}}));
      mc->add(feasst::MakeTrialGrow({
        {{"particle_type", feasst::str(particle_type)}, {"weight", "0.5"},
         {"bond", "true"}, {"mobile_site", "1"}, {"anchor_site", "0"},
         {"reference_index", ref}, {"num_steps", num}}}));
      mc->add(feasst::MakeTrialGrow({
        {{"particle_type", feasst::str(particle_type)}, {"weight", "0.5"},
         {"bond", "true"}, {"mobile_site", "0"}, {"anchor_site", "1"},
         {"reference_index", ref}, {"num_steps", num}}}));
    } else {
      mc->add(feasst::MakeTrialRotate({{"particle_type", feasst::str(particle_type)}, {"weight", "1."},
        {"tunable_param", "1."}, {"reference_index", ref}, {"num_steps", num}}));
      mc->add(feasst::MakeTrialTransfer({{"particle_type", feasst::str(particle_type)},
        {"weight", "4"}, {"reference_index", ref}, {"num_steps", num}}));
    }
  }
  mc->add(feasst::MakeCheckEnergy({{"steps_per", steps_per}, {"tolerance", "0.0001"}}));
  mc->add(feasst::MakeTuner({{"steps_per", steps_per}, {"stop_after_phase", "0"}}));
  mc->add(feasst::MakeLogAndMovie({{"steps_per", steps_per},
                                   {"file_name", "clones" + feasst::str(thread)},
                                   {"file_name_append_phase", "True"}}));
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
                                 {"num_hours_terminate", feasst::str(0.99*args.get_int("--num_procs")*args.get_double("--num_hours"))}}));
  mc->add(feasst::MakeCheckRigidBonds({{"steps_per", steps_per}}));
  return mc;
}

int main(int argc, char ** argv) {
  std::cout << feasst::version() << std::endl;
  std::cout << "args: " << args.parse(argc, argv) << std::endl;
  const int max_particles = args.get_int("--max_particles");
  const int min_particles = args.get_int("--min_particles");
  const int num_procs = args.get_int("--num_procs");
  auto windows = feasst::WindowExponential({
    {"alpha", "1.75"},
    {"num", feasst::str(num_procs)},
    {"minimum", feasst::str(min_particles)},
    {"maximum", feasst::str(max_particles)},
    {"extra_overlap", "2"}}).boundaries();
  std::cout << feasst::feasst_str(windows) << std::endl;
  auto clones = feasst::MakeClones();
  if (args.get_int("--task") == 0) {
    for (int proc = 0; proc < static_cast<int>(windows.size()); ++proc) {
      clones->add(mc(proc, windows[proc][0], windows[proc][1]));
    }
    clones->set(feasst::MakeCheckpoint({{"file_name", "checkpoint.fst"}}));
  } else {
    clones = feasst::MakeClones("checkpoint", num_procs);
  }
  clones->initialize_and_run_until_complete({{"ln_prob_file", "ln_prob.txt"}});
  std::cout << feasst::feasst_str(clones->ln_prob().values()) << std::endl;
  std::ofstream file("clones.fst");
  file << clones->serialize();
  file.close();
}
