#include "monte_carlo/include/trial_select_particle.h"

namespace feasst {

TrialSelectParticle::TrialSelectParticle(const argtype& args)
  : TrialSelect(args) {
  class_name_ = "TrialSelectParticle";
  Arguments args_(args);
  args_.dont_check();
  load_coordinates_ = args_.key("load_coordinates").dflt("true").boolean();

  // parse site
  site_ = args_.key("site").dflt("-1").integer();
  if (site_ != -1) {
    mobile_.clear();
    mobile_.add_site(0, site_);
    site_vec_ =  {site_};
  }

  set_ghost(args_.key("ghost").dflt("false").boolean());
}

class MapTrialSelectParticle {
 public:
  MapTrialSelectParticle() {
    auto obj = MakeTrialSelectParticle();
    obj->deserialize_map()["TrialSelectParticle"] = obj;
  }
};

static MapTrialSelectParticle mapper_ = MapTrialSelectParticle();

int TrialSelectParticle::random_particle(const Configuration& config,
    SelectPosition * select,
    Random * random) {
  ASSERT(group_index() >= 0, "error");
  const int num = config.num_particles(group_index());
  if (num > 0) {
    const int index = random->uniform(0, num - 1);
    const SelectGroup& ran = config.group_select(group_index());
    DEBUG("index " << group_index() << " " << index);
    DEBUG("num " << ran.num_particles());
    bool fast;
    if (site_ == - 1) {
      fast = select->replace_indices(ran.particle_index(index),
                                     ran.site_indices(index));
    } else {
      fast = select->replace_indices(ran.particle_index(index),
                                     site_vec_);
    }
    if (load_coordinates()) {
      if (!fast) select->resize();
      select->load_positions(config.particles());
    }
  } else {
    select->clear();
  }
  return num;
}

void TrialSelectParticle::ghost_particle(Configuration * config,
  SelectPosition * select) {
  ASSERT(static_cast<int>(config->ghosts().size()) > particle_type(),
    "type not recognized");
  // if no ghosts, create one
  DEBUG("particle_type: " << particle_type());
  DEBUG("nump " << config->num_particles());
  DEBUG("num ghosts " << config->ghosts()[particle_type()].num_particles());
  if (config->ghosts()[particle_type()].num_particles() == 0) {
    config->add_particle_of_type(particle_type());
    Select add;
    DEBUG("newest particle " << config->newest_particle_index());
    add.add_particle(config->newest_particle(), config->newest_particle_index());
    DEBUG("add sel: " << add.str());
    config->remove_particles(add);
    const int num_ghosts = config->ghosts()[particle_type()].num_particles();
    ASSERT(num_ghosts == 1,
      "ghost wasn't added as expected, num: " << num_ghosts);
  }
  const Select& ghost = config->ghosts()[particle_type()];
  bool fast;
  // replace indices with the last ghost as may be optimal method available
  // to delete.
  if (site_ == -1) {
    fast = select->replace_indices(ghost.particle_indices().back(),
                                   ghost.site_indices().back());
  } else {
    fast = select->replace_indices(ghost.particle_indices().back(),
                                   site_vec_);
    config->set_selection_physical(ghost, false);
    config->set_selection_physical(*select, true);
  }
  if (load_coordinates()) {
    if (!fast) {
      select->resize();
    }
    select->load_positions(config->particles());
  }
}

bool TrialSelectParticle::select(const Select& perturbed,
                                 System* system,
                                 Random * random) {
  if (is_ghost()) {
    ghost_particle(system->get_configuration(), &mobile_);
    set_probability(1.);
  } else {
    const int num = random_particle(system->configuration(), &mobile_, random);
    if (num <= 0) return false;
    set_probability(1./static_cast<double>(num));
  }
  mobile_.remove_unphysical_sites(system->configuration());
  mobile_original_ = mobile_;
  return true;
}

std::shared_ptr<TrialSelect> TrialSelectParticle::create(std::istream& istr) const {
  return std::make_shared<TrialSelectParticle>(istr);
}

TrialSelectParticle::TrialSelectParticle(std::istream& istr)
  : TrialSelect(istr) {
  // ASSERT(class_name_ == "TrialSelectParticle", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(760 == version, "mismatch version: " << version);
  feasst_deserialize(&load_coordinates_, istr);
  feasst_deserialize(&site_, istr);
  feasst_deserialize(&site_vec_, istr);
}

void TrialSelectParticle::serialize_trial_select_particle_(std::ostream& ostr) const {
  serialize_trial_select_(ostr);
  feasst_serialize_version(760, ostr);
  feasst_serialize(load_coordinates_, ostr);
  feasst_serialize(site_, ostr);
  feasst_serialize(site_vec_, ostr);
}

void TrialSelectParticle::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_select_particle_(ostr);
}

}  // namespace feasst
