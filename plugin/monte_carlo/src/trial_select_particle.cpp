#include "utils/include/serialize.h"
#include "utils/include/utils.h"  // find_in_list
#include "math/include/random.h"
#include "monte_carlo/include/trial_select_particle.h"

namespace feasst {

TrialSelectParticle::TrialSelectParticle(argtype args)
  : TrialSelectParticle(&args) {
  FEASST_CHECK_ALL_USED(args);
}

TrialSelectParticle::TrialSelectParticle(argtype * args) : TrialSelect(args) {
  class_name_ = "TrialSelectParticle";
  load_coordinates_ = boolean("load_coordinates", args, true);

  // parse site
  site_ = integer("site", args, -1);
  if (site_ != -1) {
    mobile_.clear();
    mobile_.add_site(0, site_);
    site_vec_ =  {site_};
  }

  set_ghost(boolean("ghost", args, false));
  exclude_perturbed_ = boolean("exclude_perturbed", args, false);
}

class MapTrialSelectParticle {
 public:
  MapTrialSelectParticle() {
    auto obj = MakeTrialSelectParticle();
    obj->deserialize_map()["TrialSelectParticle"] = obj;
  }
};

static MapTrialSelectParticle mapper_ = MapTrialSelectParticle();

int TrialSelectParticle::num_excluded_(const Configuration& config,
    const Select * exclude) {
  int num_excluded = 0;
  if (exclude) {
    ASSERT(is_particle_type_set(),
      "exclusion requires group to be defined by particle type to "
        << "accurately compute number of choices.");
    for (int ipart = 0; ipart < exclude->num_particles(); ++ipart) {
      DEBUG("ipart " << ipart);
      DEBUG(" particle_type() " <<  particle_type());
      if ( (particle_type() == -1 ) ||
           (config.select_particle(exclude->particle_index(ipart)).type() ==
            particle_type()) ) {
        ++num_excluded;
      }
    }
  }
  DEBUG("num_excluded " << num_excluded);
  return num_excluded;
}

int TrialSelectParticle::random_particle(const Configuration& config,
    const Select * exclude,
    Select * select,
    Random * random) {
  ASSERT(group_index() >= 0, "error");
  const int num_excluded = num_excluded_(config, exclude);
  DEBUG("grp " << group_index());
  const int num = config.num_particles(group_index());
  if (num - num_excluded > 0) {
    int attempts = 0;
    const Select * ran;
    int sel_index = -1;
    while (sel_index == -1) {
      const int index = random->uniform(0, num - 1);
      ran = const_cast<Select*>(&config.group_select(group_index()));
      DEBUG("index " << group_index() << " " << index);
      DEBUG("num " << ran->num_particles());
      if (exclude) {
        DEBUG("exclude: " << exclude->str());
        DEBUG("ran->particle_index(index) " << ran->particle_index(index));
        if (!find_in_list(ran->particle_index(index),
                          exclude->particle_indices())) {
          sel_index = index;
        }
      } else {
        DEBUG("not excluded");
        sel_index = index;
      }
      ASSERT(attempts < 1000*(num - num_excluded), "attempts:" << attempts
        << " infinite loop?");
      ++attempts;
    }
    bool fast;
    if (site_ == - 1) {
      fast = select->replace_indices(ran->particle_index(sel_index),
                                     ran->site_indices(sel_index));
    } else {
      fast = select->replace_indices(ran->particle_index(sel_index),
                                     site_vec_);
    }
    if (load_coordinates()) {
      if (!fast) select->resize_positions();
      select->load_positions(config.particles());
    }
  } else {
    select->clear();
  }
  return num - num_excluded;
}

void TrialSelectParticle::ghost_particle(Configuration * config,
  const Select * exclude,
  Select * select) {
  ASSERT(static_cast<int>(config->ghosts().size()) > particle_type(),
    "type not recognized");
  // if no ghosts, create one
  DEBUG("particle_type: " << particle_type());
  DEBUG("nump " << config->num_particles());
  DEBUG("num ghosts " << config->ghosts()[particle_type()].num_particles());
  const Select& ghosts = config->ghosts()[particle_type()];
  const int num_excluded = num_excluded_(
    const_cast<const Configuration&>(*config), exclude);
  int pindex = -1;
  if (ghosts.num_particles() <= num_excluded) {
    config->add_non_ghost_particle_of_type(particle_type());
    Select add;
    DEBUG("newest particle " << config->newest_particle_index());
    add.add_particle(config->newest_particle(), config->newest_particle_index());
    DEBUG("add sel: " << add.str());
    config->remove_particles(add);
    const int num_ghosts = ghosts.num_particles();
    ASSERT(num_ghosts > num_excluded,
      "ghost wasn't added as expected, num: " << num_ghosts);
    pindex = ghosts.num_particles() - 1;
  } else if (num_excluded > 0) {
    DEBUG("find an existing ghost that isn't excluded");
    for (int si = ghosts.num_particles() - 1; si >= 0; --si) {
      if (!find_in_list(ghosts.particle_index(si),
                        exclude->particle_indices())) {
        pindex = si;
      }
    }
  } else {
    pindex = ghosts.num_particles() - 1;
  }
  bool fast;
  // replace indices with the last ghost as may be optimal method available
  // to delete.
  DEBUG("site_ " << site_);
  DEBUG("ghosts " << ghosts.str());
  if (site_ == -1) {
    fast = select->replace_indices(ghosts.particle_index(pindex),
                                   ghosts.site_indices(pindex));
  } else {
    fast = select->replace_indices(ghosts.particle_index(pindex),
                                   site_vec_);
  }
  //config->set_selection_physical(ghosts, false);
  //config->set_selection_physical(*select, true);
  if (load_coordinates()) {
    DEBUG("fast " << fast);
    if (!fast) select->resize_positions();
    select->load_positions(config->particles());
  }
}

bool TrialSelectParticle::select(const Select& perturbed,
                                 System* system,
                                 Random * random) {
  DEBUG("is_ghost " << is_ghost());
  if (is_ghost()) {
    if (exclude_perturbed_) {
      ghost_particle(system->get_configuration(),
        const_cast<Select*>(&perturbed), &mobile_);
    } else {
      ghost_particle(system->get_configuration(), &mobile_);
    }
    set_probability_(1.);
  } else {
    int num = -1;
    if (exclude_perturbed_) {
      num = random_particle(system->configuration(),
        const_cast<Select*>(&perturbed), &mobile_, random);
    } else {
      num = random_particle(system->configuration(), &mobile_, random);
    }
    DEBUG("num " << num);
    if (num <= 0) return false;
    set_probability_(1./static_cast<double>(num));
  }
  DEBUG("selected " << mobile_.str());
  remove_unphysical_sites(system->configuration());
  ASSERT(mobile_.num_particles() > 0, "all sites shouldn't be unphysical");
  mobile_original_ = mobile_;
  DEBUG("selected " << mobile_.str());
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
  feasst_deserialize(&exclude_perturbed_, istr);
}

void TrialSelectParticle::serialize_trial_select_particle_(std::ostream& ostr) const {
  serialize_trial_select_(ostr);
  feasst_serialize_version(760, ostr);
  feasst_serialize(load_coordinates_, ostr);
  feasst_serialize(site_, ostr);
  feasst_serialize(site_vec_, ostr);
  feasst_serialize(exclude_perturbed_, ostr);
}

void TrialSelectParticle::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_select_particle_(ostr);
}

void TrialSelectParticle::select_particle(const int index,
    const Configuration& config) {
  mobile_.clear();
  const Select& select = config.group_selects()[group_index()];
  ASSERT(index < select.num_particles(), "error");
  mobile_.add_particle(select.particle_index(index), select.site_indices(index));
  mobile_.resize_positions();
  mobile_.load_positions(config.particles());
}

}  // namespace feasst
