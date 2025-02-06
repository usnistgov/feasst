#include <cmath>
#include "utils/include/arguments.h"
#include "utils/include/serialize.h"
#include "utils/include/arguments.h"
#include "math/include/accumulator.h"
#include "configuration/include/domain.h"
#include "configuration/include/configuration.h"
#include "system/include/system.h"
#include "system/include/thermo_params.h"
#include "monte_carlo/include/criteria.h"
#include "monte_carlo/include/trial_volume.h"
#include "monte_carlo/include/trial_remove.h"
#include "monte_carlo/include/trial_add.h"
#include "monte_carlo/include/monte_carlo.h"
#include "gibbs/include/gibbs_initialize.h"

namespace feasst {

FEASST_MAPPER(GibbsInitialize,);

GibbsInitialize::~GibbsInitialize() {}

void GibbsInitialize::serialize(std::ostream& ostr) const {
  Stepper::serialize(ostr);
  feasst_serialize_version(2068, ostr);
  feasst_serialize(frac_part_low_dens_, ostr);
  feasst_serialize(low_dens_config_, ostr);
  feasst_serialize(high_dens_config_, ostr);
  feasst_serialize(particle_type_, ostr);
  feasst_serialize(updates_since_adjust_, ostr);
  feasst_serialize(updates_density_equil_, ostr);
  feasst_serialize(updates_per_adjust_, ostr);
  feasst_serialize(frac_part_low_dens_tol_, ostr);
}

GibbsInitialize::GibbsInitialize(std::istream& istr) : Modify(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(2068 == version, version);
  feasst_deserialize(&frac_part_low_dens_, istr);
  feasst_deserialize(&low_dens_config_, istr);
  feasst_deserialize(&high_dens_config_, istr);
  feasst_deserialize(&particle_type_, istr);
  feasst_deserialize(&updates_since_adjust_, istr);
  feasst_deserialize(&updates_density_equil_, istr);
  feasst_deserialize(&updates_per_adjust_, istr);
  feasst_deserialize(&frac_part_low_dens_tol_, istr);
}

GibbsInitialize::GibbsInitialize(argtype * args) : Modify(args) {
  frac_part_low_dens_ = dble("fraction_particles_low_density", args, 0.15);
  particle_type_ = integer("particle_type", args, 0);
  updates_density_equil_ = integer("updates_density_equil", args, 1e7);
  updates_per_adjust_ = integer("updates_per_adjust", args, 2e7);
  frac_part_low_dens_tol_ = dble("frac_part_low_dens_tol", args, 0.05);
}
GibbsInitialize::GibbsInitialize(argtype args) : GibbsInitialize(&args) {
  feasst_check_all_used(args);
}

int GibbsInitialize::dens_config_index_(const bool high,
    const System& system) const {
  const Configuration& c0 = system.configuration(0);
  const Configuration& c1 = system.configuration(1);
  const int num0 = c0.num_particles_of_type(particle_type_);
  const int num1 = c1.num_particles_of_type(particle_type_);
  const double rho0 = static_cast<double>(num0)/c0.domain().volume();
  const double rho1 = static_cast<double>(num1)/c1.domain().volume();
  int index = 0;
  if (rho1 > rho0) {
    if (high) {
      index = 1;
    } else {
      index = 0;
    }
  } else {
    if (high) {
      index = 0;
    } else {
      index = 1;
    }
  }
  DEBUG("rho0:" << rho0 << " rho1:" << rho1 << " high dens index " << index);
  return index;
}

void GibbsInitialize::initialize(MonteCarlo * mc) {
  const System& system = mc->system();
  DEBUG("initializing");
  ASSERT(system.num_configurations() == 2, "Assumes 2 domains, but there are: "
    << system.num_configurations());
  high_dens_config_ = dens_config_index_(true, system);
  low_dens_config_ = dens_config_index_(false, system);
  DEBUG("done initializing");
}

void GibbsInitialize::update(MonteCarlo * mc) {
  DEBUG("updating");
  System * system = mc->get_system();
  Criteria * criteria = mc->get_criteria();
  ++updates_since_adjust_;
  DEBUG("serialize and input these arguments");
  const Configuration& lowd_conf = system->configuration(low_dens_config_);
  const Configuration& highd_conf = system->configuration(high_dens_config_);
  const double vol_low = lowd_conf.domain().volume();
  const double vol_high = highd_conf.domain().volume();
  const int num_low = lowd_conf.num_particles();
  const int num_high = highd_conf.num_particles();
  DEBUG("obtaining ensemble averages of the densities");
  if (!low_dens_) {
    low_dens_ = std::make_unique<Accumulator>();
  }
  if (!high_dens_) {
    high_dens_ = std::make_unique<Accumulator>();
  }
  if (updates_since_adjust_ >= updates_density_equil_) {
    low_dens_->accumulate(static_cast<double>(num_low)/vol_low);
    high_dens_->accumulate(static_cast<double>(num_high)/vol_high);
  }
  DEBUG("updates_since_adjust_ " << updates_since_adjust_);
  DEBUG("updates_per_adjust_ " << updates_per_adjust_);
  if (updates_since_adjust_ >= updates_per_adjust_) {
    INFO("Begin adjustment process. Checking if tolerances reached.");
    INFO("Estimated low and high densities: " << low_dens_->average() << " "
      << high_dens_->average());
    const double frac_low = static_cast<double>(num_low)/(num_low + num_high);
    const bool frac_tol = std::abs(frac_low - frac_part_low_dens_) <= frac_part_low_dens_tol_;
    INFO("May adjust volume to reach target fraction N in low dens:"
      << frac_part_low_dens_ << " currently: " << frac_low << " tolerance: "
      << frac_part_low_dens_tol_ << " within? " << frac_tol);
    if (frac_tol) {
      criteria->set_complete();
    } else {
      ASSERT(low_dens_->average()*high_dens_->average() != 0,
        "GibbsInitialize requires non-zero density averages.");
      const double v_low_target = frac_part_low_dens_*(num_low + num_high)/low_dens_->average();
      const double v_high_target = (1. - frac_part_low_dens_)*(num_low + num_high)/high_dens_->average();
      const double dv = v_low_target + v_high_target - vol_low - vol_high;
      INFO("Volume change:" << dv);
      system->change_volume(dv, {{"configuration", str(low_dens_config_)}});
      criteria->initialize(system);
      INFO("Energy of low density:" << criteria->current_energy(low_dens_config_));
      // HWH consider a series of trials to reach the target volume
    }

    // reset stats
    low_dens_ = std::make_unique<Accumulator>();
    high_dens_ = std::make_unique<Accumulator>();
    updates_since_adjust_ = 0;
  }
}

std::string GibbsInitialize::header(const MonteCarlo& mc) const {
  std::stringstream ss;
  return ss.str();
}

std::string GibbsInitialize::write(MonteCarlo * mc) {
  std::stringstream ss;
  return ss.str();
}

}  // namespace feasst
