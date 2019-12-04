#include "steppers/include/mean_squared_displacement.h"

namespace feasst {

class MapMeanSquaredDisplacement {
 public:
  MapMeanSquaredDisplacement() {
    MeanSquaredDisplacement().deserialize_map()["MeanSquaredDisplacement"] = MakeMeanSquaredDisplacement();
  }
};

static MapMeanSquaredDisplacement mapper_energy_check_ = MapMeanSquaredDisplacement();

MeanSquaredDisplacement::MeanSquaredDisplacement(const argtype &args)
  : Modify(args) {
  args_.init(args);
  updates_per_origin_ = args_.key("updates_per_origin").dflt("1000").integer();
  group_index_ = args_.key("group_index").dflt("0").integer();
}

void MeanSquaredDisplacement::initialize(Criteria * criteria,
    System * system,
    TrialFactory * trial_factory) {
  updates_since_origin_ = updates_per_origin_;
  system->get_configuration()->init_wrap(false);
}

void MeanSquaredDisplacement::update(Criteria * criteria,
    System * system,
    TrialFactory * trial_factory) {

  // check for new origins
  if (updates_since_origin_ == updates_per_origin_) {
    const SelectGroup& new_origin = system->configuration().group_select(group_index_);
    if (origins_.size() != 0) {
      ASSERT(new_origin.num_particles() == origins_.back().num_particles(),
        "The previous number of particles: " << origins_.back().num_particles() <<
        " does not match current number: " <<
        new_origin.num_particles() <<
        " This implementation assumes constant number of particles");
    }
    SelectPosition new_origin_coord =
      SelectPosition(new_origin, system->configuration().particles());
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
                system->configuration(),
                origins_[orig_index]);
  }

  DEBUG("updates_since_origin_" << updates_since_origin_);
  DEBUG("frames_ " << num_frames_() << " origins " << origins_.size());
}

std::string MeanSquaredDisplacement::write(Criteria * criteria,
    System * system,
    TrialFactory * trial_factory) {
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
  : Modify(istr) {
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
    const SelectPosition& old_parts) {
  if (updates >= static_cast<int>(msd_.size())) {
    msd_.resize(updates + 1);
  }
  const SelectGroup& select = config.group_select(group_index_);
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
        rsq += pow(rn - ro, 2);
      }
      msd_[updates].accumulate(rsq);
    }
  }
}

}  // namespace feasst
