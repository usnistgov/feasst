#include <fstream>
#include "utils/include/arguments.h"
#include "utils/include/serialize.h"
#include "utils/include/debug.h"
#include "configuration/include/site.h"
#include "configuration/include/particle.h"
#include "configuration/include/configuration.h"
#include "configuration/include/select.h"
#include "system/include/system.h"
#include "monte_carlo/include/trial_select_particle.h"
#include "monte_carlo/include/perturb_translate.h"
#include "monte_carlo/include/monte_carlo.h"
#include "actions/include/energy_grid.h"

namespace feasst {

EnergyGrid::EnergyGrid(argtype * args) {
  class_name_ = "EnergyGrid";
  config_ = str("config", args, "0");
  for (const std::string& rs : split(str("resolution", args, ""), ',')) {
    std::vector<double> rlist;
    for (const std::string& rs2 : split(rs, ':')) {
      rlist.push_back(str_to_double(rs2));
    }
    res_.push_back(rlist);
  }
  output_file_ = str("output_file", args, "");
}
EnergyGrid::EnergyGrid(argtype args) : EnergyGrid(&args) {
  feasst_check_all_used(args);
}
EnergyGrid::~EnergyGrid() {}
FEASST_MAPPER(EnergyGrid,);

EnergyGrid::EnergyGrid(std::istream& istr) : Action(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 2648, "mismatch version: " << version);
  feasst_deserialize(&config_, istr);
  feasst_deserialize(&output_file_, istr);
  feasst_deserialize(&res_, istr);
}

void EnergyGrid::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_action_(ostr);
  feasst_serialize_version(2648, ostr);
  feasst_serialize(config_, ostr);
  feasst_serialize(output_file_, ostr);
  feasst_serialize(res_, ostr);
}

void EnergyGrid::run(MonteCarlo * mc) {
  if (static_cast<int>(res_.size()) != 2 || output_file_.empty()) {
    WARN("No Action taken");
  }
  System * sys = mc->get_system();
  const int config_index = sys->configuration_index(config_);
  //Configuration * config = sys->get_configuration(config_index);
  const Configuration& config = sys->configuration(config_index);
  ASSERT(config.num_particles() == 2, "EnergyGrid assumes just two particles, "
    << "but there are " << config.num_particles());

  TrialSelectParticle sel;
  sel.precompute(sys);
  sel.select_particle(1, config);
  sel.set_mobile_original(sys);

  PerturbTranslate trans;
  trans.precompute(&sel, sys);
  trans.set_revert_possible(true, &sel);
  //Select select(1, config->particle(1));
  Position traj(argtype({{"dimension", str(config.dimension())}}));

  DEBUG(config.num_particles_of_type(0));
  DEBUG(config.num_particles_of_type(1));
  DEBUG(config.particle(0).site(0).position().str());
  DEBUG(config.particle(0).site(0).euler().str());
  DEBUG(config.particle(1).site(0).position().str());

  std::ofstream file(output_file_);
  file << "x,y,e" << std::endl;
  for (double x = res_[0][0]; x < res_[0][1] + res_[0][2]/2.; x += res_[0][2]) {
    traj.set_coord(0, x);
    for (double y = res_[1][0]; y < res_[1][1] + res_[1][2]/2.; y += res_[1][2]) {
      traj.set_coord(1, y);
      DEBUG("traj:" << traj.str());
      sel.select_particle(1, config);
      DEBUG("mobile pos after sel " << sel.mobile().site_positions()[0][0].str());
      trans.move(traj, sys, &sel);
      DEBUG("mobile:" << sel.mobile().str() << " pos:" << sel.mobile().site_positions()[0][0].str());
      //DEBUG("mobile:" << select.str() << " pos:" << select.site_positions()[0][0].str());
      file << x << "," << y << "," << sys->perturbed_energy(sel.mobile(), config_index) << std::endl;
      //file << x << "," << y << "," << sys->perturbed_energy(select, config_index) << std::endl;
      //file << x << "," << y << "," << sys->energy(config_index) << std::endl;
      trans.revert(sys);
      DEBUG("after revert " << config.particle(1).site(0).position().str());
    }
  }
}

}  // namespace feasst
