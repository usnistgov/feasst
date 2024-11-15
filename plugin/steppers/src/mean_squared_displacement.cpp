#include <cmath>
#include "utils/include/arguments.h"
#include "utils/include/serialize.h"
#include "math/include/accumulator.h"
#include "configuration/include/select.h"
#include "configuration/include/particle.h"
#include "configuration/include/configuration.h"
#include "system/include/system.h"
#include "monte_carlo/include/monte_carlo.h"
#include "steppers/include/mean_squared_displacement.h"

namespace feasst {

FEASST_MAPPER(MeanSquaredDisplacement,);

MeanSquaredDisplacement::MeanSquaredDisplacement(argtype * args)
  : Analyze(args) {
  updates_per_origin_ = integer("updates_per_origin", args, 1000);
  group_index_ = integer("group_index", args, 0);
}
MeanSquaredDisplacement::MeanSquaredDisplacement(argtype args)
  : MeanSquaredDisplacement(&args) {
  feasst_check_all_used(args);
}

void MeanSquaredDisplacement::initialize(MonteCarlo * mc) {
  updates_since_origin_ = updates_per_origin_;
  mc->get_system()->get_configuration(configuration_index())->init_wrap(false);
}

void MeanSquaredDisplacement::update(const MonteCarlo& mc) {
  const System& system = mc.system();
  // check for new origins
  if (updates_since_origin_ >= updates_per_origin_) {
    const Select& new_origin = configuration(system).group_select(group_index_);
    DEBUG("num particles: " << new_origin.num_particles());
    DEBUG("num sites: " << new_origin.num_sites());
    if (origins_.size() != 0) {
      ASSERT(new_origin.num_particles() == origins_.back().num_particles(),
        "The previous number of particles: " << origins_.back().num_particles() <<
        " does not match current number: " <<
        new_origin.num_particles() <<
        " This implementation assumes constant number of particles");
    }
    Select new_origin_coord =
      Select(new_origin, system.configuration().particles());
    origins_.push_back(new_origin_coord);
    updates_since_origin_ = 0;
  }
  ++updates_since_origin_;

  // add new frame for each origin
  for (int orig_index = 0;
       orig_index < static_cast<int>(origins_.size());
       ++orig_index) {
    const int start = orig_index*updates_per_origin_;
    const int elapsed = num_frames_() - start;
    update_msd_(elapsed,
                system.configuration(),
                origins_[orig_index]);
  }

  DEBUG("updates_since_origin_" << updates_since_origin_);
  DEBUG("frames_ " << num_frames_() << " origins " << origins_.size());
}

std::string MeanSquaredDisplacement::write(const MonteCarlo& mc) {
  std::stringstream ss;
  for (int frame = 0; frame < num_frames_(); ++frame) {
    ss << frame << " "
       << msd_[frame].average() << " "
       << msd_[frame].num_values()
       << std::endl;
  }
  return ss.str();
}

void MeanSquaredDisplacement::serialize(std::ostream& ostr) const {
  Stepper::serialize(ostr);
  feasst_serialize_version(732, ostr);
  feasst_serialize(updates_since_origin_, ostr);
  feasst_serialize(updates_per_origin_, ostr);
  feasst_serialize(group_index_, ostr);
  feasst_serialize_fstobj(origins_, ostr);
  feasst_serialize_fstobj(msd_, ostr);
}

MeanSquaredDisplacement::MeanSquaredDisplacement(std::istream& istr)
  : Analyze(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(732 == version, "version mismatch:" << version);
  feasst_deserialize(&updates_since_origin_, istr);
  feasst_deserialize(&updates_per_origin_, istr);
  feasst_deserialize(&group_index_, istr);
  feasst_deserialize_fstobj(&origins_, istr);
  feasst_deserialize_fstobj(&msd_, istr);
}

void MeanSquaredDisplacement::update_msd_(const int updates,
    const Configuration& config,
    const Select& old_parts) {
  if (updates >= static_cast<int>(msd_.size())) {
    msd_.resize(updates + 1);
  }
  const Select& select = config.group_select(group_index_);
  DEBUG("num: " << config.num_particles());
  for (int select_index = 0;
       select_index < select.num_particles();
       ++select_index) {
    DEBUG("sel: " << select_index);
    const int particle_index = select.particle_index(select_index);
    const Particle& part = config.select_particle(particle_index);
    for (int sel_site = 0;
         sel_site < select.num_sites(select_index);
         ++sel_site) {
      const int site_index = select.site_index(particle_index, sel_site);
      double rsq = 0;
      for (int dim = 0; dim < config.dimension(); ++dim) {
        const double ro = old_parts.site_positions()[select_index][sel_site].coord(dim);
        DEBUG(part.site(site_index).position().dimension());
        DEBUG("dim: " << dim);
        const double rn = part.site(site_index).position(dim);
        rsq += std::pow(rn - ro, 2);
      }
      msd_[updates].accumulate(rsq);
    }
  }
}

}  // namespace feasst
