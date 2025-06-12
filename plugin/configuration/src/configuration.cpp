#include "utils/include/arguments.h"
#include "utils/include/utils.h"
#include "utils/include/debug.h"
#include "utils/include/io.h"
#include "utils/include/serialize.h"
#include "math/include/constants.h"
#include "math/include/table.h"
#include "configuration/include/particle_factory.h"
#include "configuration/include/neighbor_criteria.h"
#include "configuration/include/model_params.h"
#include "configuration/include/group.h"
#include "configuration/include/file_xyz.h"
#include "configuration/include/domain.h"
#include "configuration/include/select.h"
#include "configuration/include/physical_constants.h"
#include "configuration/include/configuration.h"

namespace feasst {

Configuration::Configuration(argtype args) : Configuration(&args) {
  feasst_check_all_used(args);
}
Configuration::Configuration(argtype * args) {
  domain_ = std::make_shared<Domain>(args);
  particle_types_ = std::make_shared<ParticleFactory>();
  particle_types_->unique_particles();
  unique_types_ = std::make_shared<ParticleFactory>();
  unique_types_->unique_types();
  particles_ = std::make_shared<ParticleFactory>();
  // reset_unique_indices_();
  add(MakeGroup());  // add empty group which represents all particles

  DEBUG("parse physical constants");
  std::stringstream ss(feasst::str("physical_constants", args, "CODATA2018"));
  set_physical_constants(MakeCODATA2014()->deserialize(ss));

  DEBUG("parse types");
  std::string start;
  // if only one particle type, allow drop the subscript
  start.assign("particle_type");
  if (used(start, *args)) {
    const std::vector<std::string> ptypes = split(feasst::str(start, args), ',');
    for (int ipt = 0; ipt < static_cast<int>(ptypes.size()); ++ipt) {
      std::string value = ptypes[ipt];
      std::string name;
      if (value.find(':') != std::string::npos) {
        const std::vector<std::string> v2 = split(value, ':');
        name = v2[0];
        value = v2[1];
      }
      if (value.substr(0, 7) == "/feasst") {
        value.replace(0, 7, FEASST_INSTALL_DIR);
      }
      add_particle_type(value, name);
    }
  } else {
    int type = num_particle_types();
    std::stringstream key;
    key << start << type;
    while (used(key.str(), *args)) {
      WARN("Deprecated Configuration::particle_type[i]->particle_type");
      add_particle_type(feasst::str(key.str(), args));
      ++type;
      ASSERT(type < 1e8, "type(" << type << ") is very high. Infinite loop?");
      key.str("");
      key << start << type;
    }
  }

  DEBUG("parse adding particles of type");
  for (int type = 0; type < num_particle_types(); ++type) {
    const std::string& pname = particle_type_to_name(type);
    const std::string key = "add_num_" + pname + "_particles";
    if (used(key, *args)) {
      const int nump = integer(key, args);
      for (int ipart = 0; ipart < nump; ++ipart) {
        add_particle_of_type(type);
      }
    }
  }
  for (int type = 0; type < num_particle_types(); ++type) {
    std::stringstream key;
    key << "add_particles_of_type" << type;
    const int num = integer(key.str(), args, 0);
    if (num > 0) {
      WARN("Deprecated Configuration::add_particles_of_type[i]->add_num_[name]_particles");
    }
    for (int i = 0; i < num; ++i) {
      add_particle_of_type(type);
    }
  }

  const std::string xyz_file = feasst::str("xyz_file", args, "");
  const std::string xyz_euler_file = feasst::str("xyz_euler_file", args, "");
  if (!xyz_file.empty() && xyz_euler_file.empty()) {
    FileXYZ().load(xyz_file, this);
  } else if (xyz_file.empty() && !xyz_euler_file.empty()) {
    MakeFileXYZ({{"euler", "true"}})->load(xyz_euler_file, this);
  } else if (!xyz_file.empty() && !xyz_euler_file.empty()) {
    FATAL("cannot read both xyz and xyz_euler files.");
  }

  // Check that the dimensions of all particle types and domain match
  check_dimensions();

  DEBUG("parse groups");
  if (used("group", *args)) {
    const std::vector<std::string> grps = split(feasst::str("group", args), ',');
    for (const std::string& grp : grps) {
      args->insert({"prepend", grp});
      add(std::make_shared<Group>(args), grp);
    }
  } else {
    start = "group";
    int index = 0;
    std::stringstream key;
    key << start << index;
    while (used(key.str(), *args)) {
      const std::string name = feasst::str(key.str(), args);
      args->insert({"prepend", name});
      auto group = std::make_shared<Group>(args);
      add(group, name);
      ++index;
      key.str("");
      key << start << index;
    }
  }

  init_wrap(boolean("wrap", args, true));

  DEBUG("parse ModelParam");
  if (args->size() != 0) {
    for (std::map<std::string, std::shared_ptr<ModelParam>>::iterator iter =
         ModelParam().deserialize_map().begin();
         iter != ModelParam().deserialize_map().end(); ++iter) {
      const std::string param = iter->first;
      if (args->size() != 0) {
        if (used(param, *args)) {
          const double value = dble(param, args);
          for (int site_type = 0; site_type < num_site_types(); ++site_type) {
            set_model_param(param, site_type, value);
          }
        }
      }
      if (args->size() != 0) {
        for (int site_type = 0; site_type < num_site_types(); ++site_type) {
          std::string param_arg = param + site_type_to_name(site_type);
          if (used(param_arg, *args)) {
            set_model_param(param, site_type, dble(param_arg, args));
          }
        }
      }
      if (args->size() != 0) {
        for (int site1 = 0; site1 < num_site_types(); ++site1) {
          for (int site2 = site1; site2 < num_site_types(); ++site2) {
            std::string param_arg = param + site_type_to_name(site1) + "_" +
                                            site_type_to_name(site2);
            if (used(param_arg, *args)) {
              set_model_param(param, site1, site2, dble(param_arg, args));
            }
          }
        }
      }
      if (args->size() != 0) {
        std::string param_arg = param + "_mixing_file";
        if (used(param_arg, *args)) {
          WARN("Deprecated argument [parameter]_mixing_file. " <<
            "Use mixing_file instead.");
          DEBUG("param " << param << " param_arg " << param_arg);
          std::vector<std::string> names;
          site_type_names(&names);
          set_model_param(param, feasst::str(param_arg, args), &names);
        }
      }
    }
  }
  const std::string model_param_file = feasst::str("model_param_file", args, "");
  if (!model_param_file.empty()) {
    std::vector<std::string> names;
    site_type_names(&names);
    set_model_param(model_param_file, &names);
  }

  // Set after param overrides
  if (boolean("set_cutoff_min_to_sigma", args, false)) {
    unique_types_->set_cutoff_min_to_sigma();
  }
  name_ = feasst::str("name", args, "0");
}

void Configuration::add_particle_type(const std::string file_name,
    std::string name) {
  if (name.empty()) {
    name = std::to_string(num_particle_types());
    DEBUG("setting default name for part:" << file_name << " as: " << name);
  }
  ASSERT(!find_in_list(name, type_to_name_), "particle type name:" << name
    << " already exists in Configuration. Particle names must be unique.");
  DEBUG("adding type: " << file_name << " of name: " << name);
  ASSERT(particles_->num() == 0, "types cannot be added after particles");
  particle_types_->add(file_name);
  unique_types_->add(file_name, name);
  ghosts_.push_back(std::make_shared<Select>());
  ASSERT(ghosts_.back()->is_group_empty(), "no ghosts in brand new type");
  type_to_file_.push_back(file_name);
  type_to_name_.push_back(name);
  num_particles_of_type_.push_back(0);
}

void Configuration::add_(const Particle particle) {
  Particle part = particle;
  particles_->add(part);
  for (std::shared_ptr<Select> select : group_selects_) {
    add_to_selection_(particles_->num() - 1, select.get());
  }
  position_tracker_(particles_->num() - 1);
}

void Configuration::add_non_ghost_particle_of_type(const int type) {
  Particle part = particle_types_->particle(type);
  part.erase_bonds();
  add_(part);
  newest_particle_index_ = particles_->num() - 1;
  ++num_particles_of_type_[type];
}

void Configuration::add_particle_of_type(const int type) {
  DEBUG("adding type: " << type);
  ASSERT(type < num_particle_types(), "type(" << type << ") is not allowed "
    << "when there are only " << num_particle_types() << " particle types");
  // decide whether to add a ghost particle or a new one
  DEBUG("numg " << ghosts_[type]->num_particles());
  if (ghosts_[type]->num_particles() == 0) {
    add_non_ghost_particle_of_type(type);
  } else {
    const int index = ghosts_[type]->particle_index(0);
    ghosts_[type]->remove_particle(index);
    for (std::shared_ptr<Select> select : group_selects_) {
      add_to_selection_(index, select.get());
    }
    newest_particle_index_ = index;
    ++num_particles_of_type_[type];
  }
}

void Configuration::remove_particle_(const int particle_index) {
  // HWH optimization: this randomly generated string is expensive.
  // reset_unique_indices_();
  const int type = particles_->particle(particle_index).type();
  --num_particles_of_type_[type];
  ghosts_[type]->add_particle(select_particle(particle_index), particle_index);
  DEBUG("type " << type);
  DEBUG("particle index " << particle_index);
  DEBUG("num particles " << num_particles());
  for (std::shared_ptr<Select> select : group_selects_) {
    select->remove_particle(particle_index);
  }
}

void Configuration::remove_particles(const Select& selection) {
  ASSERT(selection.num_particles() > 0, "no selection");
  for (int index = selection.num_particles()  - 1;
       index >= 0;
       --index) {
    remove_particle_(selection.particle_index(index));
  }
}

void Configuration::remove_particle(const Select& selection) {
  ASSERT(selection.num_particles() == 1, "assumes 1 particle selected");
  remove_particles(selection);
}

void Configuration::displace_particle_(const int particle_index,
                                       const Position &displacement) {
  particles_->displace(particle_index, displacement);
  position_tracker_(particle_index);
}

void Configuration::displace_site_(const int particle_index,
                                   const int site_index,
                                   const Position &displacement) {
  Position pos = displacement;
  pos.add(select_particle(particle_index).site(site_index).position());
  particles_->replace_position(particle_index, site_index, pos);
  position_tracker_(particle_index, site_index);
}

void Configuration::displace_particles(const Select& selection,
                                       const Position &displacement) {
  ASSERT(selection.num_particles() > 0, "no selection");
  for (int particle_index : selection.particle_indices()) {
    displace_particle_(particle_index, displacement);
  }
}

void Configuration::displace_particle(const Select& selection,
                                      const Position &displacement) {
  ASSERT(selection.num_particles() == 1, "assumes 1 particle selected");
  displace_particles(selection, displacement);
}

void Configuration::replace_position(const Select& select,
                                     const Particle& replacement) {
  ASSERT(select.num_particles() == 1, "only one particle can be selected");
  replace_position_(select.particle_index(0), replacement);
}

void Configuration::replace_position_(const int particle_index,
                                      const Particle& replacement) {
  ASSERT(select_particle(particle_index).num_sites() ==
         replacement.num_sites(), "size error");
  particles_->replace_position(particle_index, replacement);
  position_tracker_(particle_index);
}

void Configuration::replace_position_(const int particle_index,
                                      const int site_index,
                                      const Position& replacement) {
  particles_->replace_position(particle_index, site_index, replacement);
  position_tracker_(particle_index, site_index);
}

void Configuration::add(std::shared_ptr<Group> group, std::string name) {
  ASSERT(group->is_empty() || particle_types_->num() != 0,
    "add groups after particle types");
  if (name == "") {
    std::stringstream ss;
    ss << num_groups();
    name.assign(ss.str());
  }
  Select group_select;
  group->name_to_index(unique_types());
  group->add_property(name, 0.);
  group_select.set_group(group);
  init_selection_(&group_select);
  group_selects_.push_back(std::make_shared<Select>(group_select));
}

void Configuration::position_tracker_(const int particle_index,
                                      const int site_index) {
  ASSERT(site_index >= 0, "index error");
  DEBUG("update selection");
  for (std::shared_ptr<Select> select : group_selects_) {
    ASSERT(!select->group().is_spatial(), "implement updating of groups");
  }
}

void Configuration::position_tracker_(const int particle_index) {
  for (int site_index = 0;
       site_index < select_particle(particle_index).num_sites();
       ++site_index) {
    position_tracker_(particle_index, site_index);
  }
}

void Configuration::position_tracker_() {
  for (int index : selection_of_all().particle_indices()) {
    position_tracker_(index);
  }
}

void Configuration::set(std::shared_ptr<Domain> domain) {
  domain_ = domain;
  position_tracker_();
}

/// HWH add check .. domain, positions, particles, etc
void Configuration::check() const {
  particles_->check();
  particle_types_->check();
  unique_types_->check();
  // selection_of_all().check();

  ASSERT(particle_types_->num() == num_particle_types(), "er");
  ASSERT(unique_types_->num_sites() == num_site_types(), "er");

  // check that the first group is all particles in the configuration.
  ASSERT(static_cast<int>(group_selects_[0]->num_particles())
    == num_particles(),
    "The number of particles in the first group(" <<
    group_selects_[0]->num_particles() << ") is not equal to the number of " <<
    "particles: " << num_particles());
  for (int index = 0; index < num_particles(); ++index) {
    ASSERT(static_cast<int>(group_selects_[0]->site_indices(index).size()) ==
      particle(index).num_sites(), "size error");
  }

  // check number of particle types
  ASSERT(
    static_cast<int>(num_particles_of_type_.size()) == num_particle_types(),
    "size error");

  // check that a particle is not simultaneously a ghost and a real particle
  for (const std::shared_ptr<Select>& ghost : ghosts_) {
    ASSERT(!group_selects_[0]->is_overlap(*ghost),
      "ghost particle cannot also be real");
  }

  // check that a particle is not listed as a ghost twice
  for (const std::shared_ptr<Select>& ghost : ghosts_) {
    ASSERT(!has_duplicate(ghost->particle_indices()),
      "the same particle cannot be listed as a ghost twice");
  }

  model_params().check();
}

bool Configuration::are_all_sites_physical() const {
  for (const Particle& part : particles_->particles()) {
    for (const Site& site : part.sites()) {
      if (!site.is_physical()) {
        return false;
      }
    }
  }
  return true;
}

// void Configuration::check_id_(const Select& select) const {
//   check_id_(select.unique_id());
// }
// void Configuration::check_id_(const std::string id) const {
//   ASSERT(id.empty() or id == unique_indices_,
//     "If selection is obtained from a configuration, it will have a unique" <<
//     "identifier attached to it. The id of the selection and the " <<
//     "configuration were fount to be inconsistent, likely due to an " <<
//     "expired selection (e.g., particles were removed after selection)");
// }

void Configuration::update_positions(
    const std::vector<std::vector<double> > coords) {
  if (coords.size() == 0) {
    return;
  }
  ASSERT(static_cast<int>(coords.size()) == num_sites(), "the number of " <<
    "coordinates provided: " << coords.size() << " does not match the number"
    << " of sites: " << num_sites());
  DEBUG("dimension: " << dimension());
  DEBUG("num sites " << coords.size());
  DEBUG("num sites in config  " << num_sites());
  DEBUG("num particles in config " << num_particles());
  ASSERT(static_cast<int>(coords[0].size()) == dimension(), "the dimensions: "
    << coords[0].size() << " of the coordinates do not match the dimensions: "
    << dimension() << " of the configuration.");
  Position position;
  int iter_site = 0;
  for (int part_index : group_selects_[0]->particle_indices()) {
    Particle part = select_particle(part_index);
    DEBUG("part_index " << part_index);
    for (int site_index = 0;
         site_index < part.num_sites();
         ++site_index) {
      position.set_vector(coords[iter_site]);
      Site site = part.site(site_index);
      site.set_position(position);
      part.set_site(site_index, site);
      ++iter_site;
    }
    replace_position_(part_index, part);
  }
}

void Configuration::update_positions(
    const std::vector<std::vector<double> > coords,
    const std::vector<std::vector<double> > eulers) {
  update_positions(coords);
  ASSERT(dimension() == 3, "Eulers require 3 dimensions.");
  Euler euler;
  int iter_site = 0;
  for (int part_index : group_selects_[0]->particle_indices()) {
    Particle * part = get_particles_()->get_particle(part_index);
    for (int site_index = 0;
         site_index < part->num_sites();
         ++site_index) {
      Site * site = part->get_site(site_index);
      euler.set(eulers[iter_site][0],
                eulers[iter_site][1],
                eulers[iter_site][2]);
      site->set_euler(euler);
      ++iter_site;
    }
  }
}

void Configuration::displace(const Select& selection,
                             const Position &displacement) {
  for (int sp = 0; sp < selection.num_particles(); ++sp) {
    const int particle_index = selection.particle_index(sp);
    for (const int& site_index : selection.site_indices(sp)) {
      displace_site_(particle_index, site_index, displacement);
    }
  }
}

void Configuration::add_to_selection_(const int particle_index,
                                      Select * select) const {
  const Particle& part = select_particle(particle_index);
  const Group& group = select->group();
  if (group.is_in(part, particle_index)) {
    select->add_particle(particle_index, group.site_indices(part));
  }
}

void Configuration::update_selection_(const int particle_index,
                                      Select * select) const {
  const Particle& part = select_particle(particle_index);
  const Group& group = select->group();
  if (group.is_in(part, particle_index)) {
    select->add_particle(particle_index, group.site_indices(part), true);
  } else {
    select->remove_particle(particle_index);
  }
}

void Configuration::init_selection_(Select * group_select) const {
  if (num_particles() > 0) {
    for (int part_index : selection_of_all().particle_indices()) {
      add_to_selection_(part_index, group_select);
    }
  }
}

void Configuration::update_positions(const Select& select,
                                     const bool no_wrap) {
  int pindex = 0;
  for (int particle_index : select.particle_indices()) {
    int sindex = 0;
    for (int site_index : select.site_indices(pindex)) {
      // DEBUG(select.site_properties()[pindex][sindex].str());
      replace_properties_(particle_index,
                          site_index,
                          select.site_properties()[pindex][sindex]);
      replace_position_(particle_index,
                        site_index,
                        select.site_positions()[pindex][sindex]);
      ++sindex;
    }

    // check if particle needs to be wrapped, for whole-particle updates only.
    if (!no_wrap) {
      if (select.num_sites(pindex) ==
          select_particle(particle_index).num_sites()) {
        DEBUG("wrapping");
        wrap_particle(particle_index);
      }
    }
    ++pindex;
  }
  if (select.is_anisotropic()) {
    int pindex = 0;
    for (int particle_index : select.particle_indices()) {
      int sindex = 0;
      for (int site_index : select.site_indices(pindex)) {
        DEBUG("particle_index " << particle_index);
        DEBUG("site_index " << site_index);
        Particle * part = particles_->get_particle(particle_index);
        Site * site = part->get_site(site_index);
        const int site_type = site->type();
        const int particle_type = part->type();
        if (unique_type(particle_type, site_type).is_anisotropic()) {
          DEBUG(select.site_eulers().size());
          DEBUG(select.site_eulers()[pindex].size());
          DEBUG("pindex " << pindex);
          DEBUG("sindex " << sindex);
          ASSERT(sindex < static_cast<int>(select.site_eulers()[pindex].size()),
            "error");
          site->set_euler(select.site_eulers()[pindex][sindex]);
        }
        ++sindex;
      }
      ++pindex;
    }
  }
}

int Configuration::particle_type_to_group(const int particle_type) const {
  int index;
  if (num_particle_types() == 1) {
    return 0;
  }
  DEBUG("group_stored.. " << feasst_str(group_store_particle_type_));
  DEBUG("num_groups " << num_groups());
  DEBUG("particle_type " << particle_type);
  if (!find_in_list(particle_type, group_store_particle_type_, &index)) {
    return -1;
  }
  return group_store_group_index_[index];
}

int Configuration::particle_type_to_group_create(const int particle_type) {
  const int grp = particle_type_to_group(particle_type);
  if (grp != -1) {
    return grp;
  }
  add(MakeGroup({{"particle_type", particle_type_to_name(particle_type)}}));
  //add(MakeGroup({{"particle_type", feasst::str(HEREITISTHEISSUEparticle_type)}}));
  group_store_particle_type_.push_back(particle_type);
  group_store_group_index_.push_back(num_groups() - 1);
  const int index = static_cast<int>(group_store_group_index_.size()) - 1;
  DEBUG("index " << index);
  DEBUG("group_stored p.. " << feasst_str(group_store_particle_type_));
  DEBUG("group_stored g.. " << feasst_str(group_store_group_index_));
  return group_store_group_index_[index];
}

int Configuration::num_particles(const int group) const {
  TRACE("pn " << particles_->num());
  TRACE("gn " << num_ghosts_());
  TRACE("group " << group);
  if (group == 0) {
    TRACE("here " << particles_->num() - num_ghosts_());
    const int num = particles_->num() - num_ghosts_();
    ASSERT(num >= 0, "error");
    return num;
  }
  TRACE("here " << group_selects_[group]->num_particles());
  const int num = group_selects_[group]->num_particles();
  ASSERT(num >= 0, "error");
  return num;
}

int Configuration::num_sites(const int group) const {
  if (group != 0) FATAL("not implemented");

  int num_ghost_sites = 0;
  for (const std::shared_ptr<Select>& ghost : ghosts_) {
    num_ghost_sites += ghost->num_sites();
  }
  DEBUG("num ghost sites " << num_ghost_sites);
  return particles_->num_sites() - num_ghost_sites;
}

int Configuration::num_ghosts_() const {
  int num = 0;
  for (const std::shared_ptr<Select>& select : ghosts_) {
    num += select->num_particles();
  }
  return num;
}

void Configuration::revive(const Select& selection) {
  for (int particle_index : selection.particle_indices()) {
    DEBUG("reviving particle_index: " << particle_index);
    const Particle& part = select_particle(particle_index);
    const int type = part.type();
    DEBUG("ghosts " << feasst_str(ghosts_[type]->particle_indices()));
    ASSERT(find_in_list(particle_index, ghosts_[type]->particle_indices()),
      "attempting to revive a particle that isn't a ghost");
    ++num_particles_of_type_[type];
    DEBUG("ghost particles " << ghosts_[type]->num_particles());
    ghosts_[type]->remove_particle(particle_index);
    for (std::shared_ptr<Select> select : group_selects_) {
      add_to_selection_(particle_index, select.get());
    }
    position_tracker_(particle_index);
  }
}

const Particle& Configuration::particle(const int index) const {
  const int particle_index = selection_of_all().particle_index(index);
  return particles_->particle(particle_index);
}

Particle Configuration::particle(const int index,
                                 const int group) const {
  const Select& select_group = *group_selects_[group];
  const int particle_index = select_group.particle_index(index);
  Particle part = particles_->particle(particle_index);
  select_group.group().remove_sites(&part);
  return part;
}

void Configuration::add(std::shared_ptr<ModelParam> param) {
  for (const Particle& particle : unique_types().particles()) {
    param->add(particle);
  }
  unique_types_->add(param);
}

int Configuration::num_particles_of_type(const int type) const {
  ASSERT(type < num_particle_types(), "type: " << type << " >= num types: "
    << num_particle_types());
  if (type == -1) {
    return num_particles();
  }
  const int num = num_particles_of_type_[type];
  ASSERT(num >= 0, "unphysical number(" << num <<")");
  ASSERT(num <= num_particles(), "unphysical number of particles of type(" <<
    num << ") when number of particles is: " << num_particles());
  return num;
}

void Configuration::wrap_particle(const int particle_index) {
  if (wrap_) {
    const Position& site0_position =
      select_particle(particle_index).site(0).position();
    const Position& pbc_shift = domain_->shift_opt(site0_position);
    DEBUG("site0_position " << site0_position.str());
    DEBUG("pbc " << pbc_shift.str());
    if (pbc_shift.squared_distance() > NEAR_ZERO) {
      displace_particle_(particle_index, pbc_shift);
      DEBUG("new pos " <<
        select_particle(particle_index).site(0).position().str());
    }
  }
}

void Configuration::set_selection_physical(const Select& select,
    const bool phys) {
  for (int sp_index = 0;
       sp_index < static_cast<int>(select.particle_indices().size());
       ++sp_index) {
    const std::vector<int>& site_indices = select.site_indices()[sp_index];
    for (int ss_index = 0;
         ss_index < static_cast<int>(site_indices.size());
         ++ss_index) {
      particles_->set_site_physical(
        select.particle_indices()[sp_index],
        site_indices[ss_index],
        phys);
    }
  }
}

bool Configuration::is_equal(const Configuration& configuration,
                             const double tolerance) const {
  // check particles/sites of non-ghosts.
  if (!selection_of_all().is_equal(configuration.selection_of_all())) {
    DEBUG("unequal selection");
    return false;
  }

  // check position of first particle.
  if (num_particles() > 0) {
    for (int pindex = 0; pindex < num_particles(); ++pindex) {
      const Particle p1 = particle(pindex);
      const Particle p2 = configuration.particle(pindex);
      for (int is = 0; is < p1.num_sites(); ++is) {
        if (!p1.site(is).position().is_equal(p2.site(is).position(),
                                             tolerance)) {
          DEBUG("unequal site" << is << " positions: "
            << p1.site(is).position().str() << " vs "
            << p2.site(is).position().str());
          return false;
        }
      }
    }
  }
  return true;
}

int Configuration::max_sites_in_any_particle() const {
  int mx = 0;
  for (int type = 0; type < num_particle_types(); ++type) {
    mx = std::max(mx, particle_type(type).num_sites());
  }
  return mx;
}

void Configuration::set_site_type(const int particle_type,
                                  const int site,
                                  const int site_type) {
  // Check if cell needs to be updated with changing type
  for (std::shared_ptr<Select> select : group_selects_) {
    if (find_in_list(site_type, select->group().site_types())) {
      ERROR("check if groups need to be updated with changing type");
    }
  }
  particle_types_->set_site_type(particle_type, site, site_type);
  for (int particle = 0; particle < particles_->num(); ++particle) {
    if (particles_->particle(particle).type() == particle_type) {
      particles_->set_site_type(particle, site, site_type);
    }
  }
}

std::vector<int> Configuration::num_sites_of_type(
    const Select& selection) const {
  std::vector<int> num;
  num_sites_of_type(selection, &num);
  return num;
}

void Configuration::num_sites_of_type(const Select& selection,
    std::vector<int> * num) const {
  if (static_cast<int>(num->size()) != num_site_types()) {
    num->resize(num_site_types());
  }
  std::fill(num->begin(), num->end(), 0);
  for (int select_index = 0;
       select_index < selection.num_particles();
       ++select_index) {
    const int part_index = selection.particle_index(select_index);
    const Particle& part = select_particle(part_index);
    for (int site_index : selection.site_indices(select_index)) {
      const Site& site = part.site(site_index);
      if (site.is_physical()) {
        (*num)[site.type()] += 1;
      }
    }
  }
}

void Configuration::set_side_lengths(const Position& sides) {
  if (domain_) {
    domain_->set_side_lengths(sides);
  } else {
    WARN("site lengths were attempted to be set without a domain");
  }
}

std::string Configuration::status_header(const std::string append) const {
  std::stringstream ss;
  ss << domain().status_header(append);
  for (int type = 0; type < num_particle_types(); ++type) {
    ss << ",num_particles_" << particle_type_to_name(type) << append;
  }
  return ss.str();
}

std::string Configuration::status() const {
  std::stringstream ss;
  ss << domain().status();
  for (int type = 0; type < num_particle_types(); ++type) {
    ss << "," << num_particles_of_type(type);
  }
  return ss.str();
}

void Configuration::serialize(std::ostream& ostr) const {
  feasst_serialize_version(7200, ostr);
  feasst_serialize(std::string(FEASST_VERSION), ostr);
  feasst_serialize(particle_types_, ostr);
  feasst_serialize(unique_types_, ostr);
  feasst_serialize(particles_, ostr);
  feasst_serialize(domain_, ostr);
  feasst_serialize(group_selects_, ostr);
  feasst_serialize(group_store_particle_type_, ostr);
  feasst_serialize(group_store_group_index_, ostr);
  feasst_serialize(ghosts_, ostr);
  feasst_serialize(type_to_file_, ostr);
  feasst_serialize(num_particles_of_type_, ostr);
  feasst_serialize(wrap_, ostr);
  feasst_serialize(num_cell_lists_, ostr);
  feasst_serialize(neighbor_criteria_, ostr);
  feasst_serialize(name_, ostr);
  feasst_serialize_endcap("Configuration", ostr);
  DEBUG("size: " << ostr.tellp());
}

Configuration::Configuration(std::istream& istr) {
  const int config_version = feasst_deserialize_version(istr);
  ASSERT(config_version >= 7199 && config_version <= 7200,
    "unrecognized config_version: " << config_version);
  std::string checkpoint_version;
  feasst_deserialize(&checkpoint_version, istr);
  if (checkpoint_version != FEASST_VERSION) {
    WARN("version of checkpoint: " << checkpoint_version << " is not the " <<
         "same as current version: " << FEASST_VERSION);
  }
//  feasst_deserialize(particle_types_, istr);
// HWH for unknown reasons, this function template does not work.
  {
    int existing;
    istr >> existing;
    if (existing != 0) {
      particle_types_ = std::make_shared<ParticleFactory>(istr);
    }
  }
  {
    int existing;
    istr >> existing;
    if (existing != 0) {
      unique_types_ = std::make_shared<ParticleFactory>(istr);
    }
  }
  {
    int existing;
    istr >> existing;
    if (existing != 0) {
      particles_ = std::make_shared<ParticleFactory>(istr);
    }
  }
  // HWH for unknown reasons, this doesn't work
  // feasst_deserialize(domain_, istr);
  { int existing;
    istr >> existing;
    if (existing != 0) domain_ = std::make_shared<Domain>(istr);
  }
  // feasst_deserialize(group_selects_, istr);
//  HWH for unknown reasons, this function template does not work.
  {
    int dim1;
    istr >> dim1;
    group_selects_.resize(dim1);
    for (int index = 0; index < dim1; ++index) {
      int existing;
      istr >> existing;
      if (existing != 0) {
        group_selects_[index] = std::make_shared<Select>(istr);
      }
    }
  }
  feasst_deserialize(&group_store_particle_type_, istr);
  feasst_deserialize(&group_store_group_index_, istr);
  // feasst_deserialize(ghosts_, istr);
//  HWH for unknown reasons, this function template does not work.
  {
    int dim1;
    istr >> dim1;
    ghosts_.resize(dim1);
    for (int index = 0; index < dim1; ++index) {
      int existing;
      istr >> existing;
      if (existing != 0) {
        ghosts_[index] = std::make_shared<Select>(istr);
      }
    }
  }
  feasst_deserialize(&type_to_file_, istr);
  feasst_deserialize(&num_particles_of_type_, istr);
  feasst_deserialize(&wrap_, istr);
  feasst_deserialize(&num_cell_lists_, istr);
//  feasst_deserialize(neighbor_criteria_, istr);
// HWH for unknown reasons, this function template does not work.
  {
    int dim1;
    istr >> dim1;
    neighbor_criteria_.resize(dim1);
    for (int index = 0; index < dim1; ++index) {
      int existing;
      istr >> existing;
      if (existing != 0) {
        neighbor_criteria_[index] = std::make_shared<NeighborCriteria>(istr);
      }
    }
  }
  if (config_version >= 7200) {
    feasst_deserialize(&name_, istr);
  }
  feasst_deserialize_endcap("Configuration", istr);
}

void Configuration::copy_particles(const Configuration& config,
    const bool add_missing) {
  for (int ipart = 0; ipart < num_particles(); ++ipart) {
    const int type1 = particle(ipart).type();
    const int nsite1 = particle(ipart).num_sites();
    const int type2 = config.particle(ipart).type();
    const int nsite2 = config.particle(ipart).num_sites();
    ASSERT(type1 == type2, "type1: " << type1 << " != type2 " << type2);
    ASSERT(nsite1 == nsite2, "nsite1: " << nsite1 << " != nsite2 " << nsite2);
  }

  if (add_missing) {
    for (int extra = num_particle_types();
         extra < config.num_particle_types();
         ++extra) {
      add_particle_type(config.type_to_file_name(extra));
    }
    for (int extra = num_particles(); extra < config.num_particles(); ++extra) {
      add_particle_of_type(config.particle(extra).type());
    }
  }
  ASSERT(num_particles() == config.num_particles(),
    "number of particles in config " << num_particles() <<
    " does not match number of given config " << config.num_particles());

  DEBUG("replace positions");
  int part2 = 0;
  for (const int part1 : selection_of_all().particle_indices()) {
    replace_position_(part1, config.particle(part2));
    ++part2;
  }
}

int Configuration::dimension() const { return domain().dimension(); }

void Configuration::synchronize_(const Configuration& config,
    const Select& perturbed) {
  DEBUG(perturbed.str());
  Select sync_sel(perturbed, *(config.particles_));
  update_positions(sync_sel, true);
//  for (int spindex = 0; spindex < perturbed.num_particles(); ++spindex) {
//    const int part_index = perturbed.particle_index(spindex);
//    one_site_select_.set_particle(0, part_index);
//    DEBUG("part_index " << part_index);
//    const Particle& part = config.select_particle(part_index);
//    for (int sindex : perturbed.site_indices(spindex)) {
////      Site * site = particles_->get_particle(part_index)->get_site(sindex);
////      *site = part.site(sindex);
////      position_tracker_(part_index, sindex);
////    }
//    *particles_->get_particle(part_index) = part;
//    position_tracker_(part_index);
//  }
}

void Configuration::set_particle_type(const int ptype,
                                      const Select& select) {
  for (int particle_index : select.particle_indices()) {
    Particle * part = particles_->get_particle(particle_index);
    --num_particles_of_type_[part->type()];
    part->set_type(ptype);
    ++num_particles_of_type_[part->type()];
    ASSERT(part->num_sites() == particle_type(ptype).num_sites(),
      "cannot morph into particle with different number of sites");
    for (int isite = 0; isite < part->num_sites(); ++isite) {
      part->get_site(isite)->set_type(particle_type(ptype).site(isite).type());
    }
    for (std::shared_ptr<Select> sel : group_selects_) {
      update_selection_(particle_index, sel.get());
    }
    // HWH doesn't update type-based cell lists, groups, etc.
  }
}

void Configuration::change_volume(const double delta_volume,
    argtype args) {
  change_volume(delta_volume, &args);
  feasst_check_all_used(args);
}

void Configuration::change_volume(const double delta_volume,
    argtype * args) {
  DEBUG("delta_volume " << delta_volume);
  ASSERT(domain().volume() + delta_volume > 0,
    "delta_volume " << delta_volume << " too large for volume "
    << domain().volume());
  ASSERT(wrap_, "positions must be wrapped before scaling for volume change");
  double factor = 1. + delta_volume/domain_->volume();
  const int dimen = integer("dimension", args, -1);
  if (dimen == -1) {
    factor = std::pow(factor, 1./dimension());
    for (int dim = 0; dim < dimension(); ++dim) {
      domain_->set_side_length(dim, domain_->side_length(dim)*factor);
    }
  } else {
    domain_->set_side_length(dimen, domain_->side_length(dimen)*factor);
  }
  if (boolean("scale_particles", args, true)) {
    particles_->scale_particle_positions(dimen, factor);
  }
  position_tracker_();
}

int Configuration::group_index(const std::string& name) const {
  for (int index = 0; index < static_cast<int>(group_selects_.size());
       ++index) {
    const Select& sel = *group_selects_[index];
    if (sel.group().has_property(name)) {
      return index;
    }
  }
  FATAL("There is no group with name: " << name);
}

std::string Configuration::str() const {
  std::stringstream ss;
  for (const Particle& part : particles_->particles()) {
    for (const Site& site : part.sites()) {
      ss << site.position().str() << std::endl;
    }
  }
  return ss.str();
}

void Configuration::check_dimensions() const {
  int dim = 0, dim_prev = 0;
  for (const Particle& part : particle_types_->particles()) {
    for (const Site& site : part.sites()) {
      dim = site.position().dimension();
      if (dim_prev != 0) {
        ASSERT(dim == dim_prev, "Particles have a mix of 2D and 3D");
      }
      dim_prev = dim;
    }
  }
  ASSERT(domain_->dimension() == 0 || domain_->dimension() == dim || dim == 0,
    "Dimension of particles: " << dim << " does not match dimension of domain: "
    << domain_->dimension());
}

int Configuration::site_type_to_particle_type(const int site_type) const {
  return unique_types().site_type_to_particle_type(site_type);
}

const Site& Configuration::unique_type(const int ptype, const int stype) const {
  return unique_types().unique_type(ptype, stype);
}

std::vector<std::vector<int> >
    Configuration::num_site_types_per_particle_type() const {
  int prev = 0;
  std::vector<std::vector<int> > nstppt(num_particle_types());
  for (int ptype = 0; ptype < num_particle_types(); ++ptype) {
    const Particle& part = particle_types().particle(ptype);
    const int num_site_types = unique_types().particle(ptype).num_sites();
    std::vector<int> st(num_site_types+prev);
    for (int stype = 0; stype < num_site_types; ++stype) {
      st[stype + prev] += part.num_sites_of_type(stype + prev);
    }
    prev += num_site_types;
    nstppt[ptype] = st;
  }
  return nstppt;
}

const PhysicalConstants& Configuration::physical_constants() const {
  return model_params().physical_constants();
}

void Configuration::add(std::shared_ptr<NeighborCriteria> neighbor_criteria) {
  neighbor_criteria_.push_back(neighbor_criteria);
  neighbor_criteria_.back()->name_to_index(*unique_types_);
}

const NeighborCriteria& Configuration::neighbor_criteria(
    const int index) const {
  return *neighbor_criteria_[index];
}

const std::vector<std::shared_ptr<NeighborCriteria> >&
  Configuration::neighbor_criteria() const { return neighbor_criteria_; }

NeighborCriteria * Configuration::get_neighbor_criteria(const int index) {
  return neighbor_criteria_[index].get();
}

int Configuration::num_particle_types() const { return particle_types_->num(); }
int Configuration::num_site_types() const { return unique_types_->num_sites(); }
int Configuration::num_bond_types() const { return unique_types_->num_bonds(); }
int Configuration::num_angle_types() const {
  return unique_types_->num_angles(); }
int Configuration::num_dihedral_types() const {
  return unique_types_->num_dihedrals(); }
const Particle& Configuration::particle_type(const int type) const {
  return particle_types_->particle(type); }
const ParticleFactory& Configuration::particle_types() const {
  return *particle_types_;
}
const ModelParams& Configuration::model_params() const {
  return unique_types_->model_params(); }

void Configuration::set_model_param(const std::string name,
                     const int site_type,
                     const double value) {
  unique_types_->set_model_param(name, site_type, value); }

void Configuration::set_model_param(const std::string name,
                     const int site_type1,
                     const int site_type2,
                     const double value) {
  unique_types_->set_model_param(name, site_type1, site_type2, value); }

void Configuration::set_model_param(const std::string name,
                                    const std::string filename,
                                    std::vector<std::string> * site_type_names) {
  unique_types_->set_model_param(name, filename, site_type_names); }

void Configuration::set_model_param(const std::string filename,
                                    std::vector<std::string> * site_type_names) {
  unique_types_->set_model_param(filename, site_type_names); }

void Configuration::add_model_param(const std::string name,
                     const double value) {
  unique_types_->add_model_param(name, value); }

void Configuration::add_or_set_model_param(const std::string name,
                            const double value) {
  unique_types_->add_or_set_model_param(name, value); }

void Configuration::set_physical_constants(
    std::shared_ptr<PhysicalConstants> constants) {
  unique_types_->set_physical_constants(constants); }

const ParticleFactory& Configuration::unique_types() const {
  return *unique_types_;
}

/// Return the unique type by individual particle.
const Particle& Configuration::unique_type(const int type) const {
  return unique_types_->particle(type); }

/// Add or set the property of a site in a particle type.
void Configuration::add_or_set_particle_type_site_property(
    const std::string name, const double value, const int particle_type,
    const int site_index) {
  particle_types_->add_or_set_site_property(name, value, particle_type,
                                            site_index);
}

const ParticleFactory& Configuration::particles() const { return *particles_; }

ParticleFactory * Configuration::get_particles_() { return particles_.get(); }

const Particle& Configuration::select_particle(const int index) const {
  return particles_->particle(index); }

void Configuration::set_property(const std::string name,
    const double value,
    const int particle_index) {
  particles_->set_property(name, value, particle_index); }

void Configuration::add_site_property(const std::string name,
    const double value,
    const int particle_index,
    const int site_index) {
  particles_->add_site_property(name, value, particle_index, site_index); }

void Configuration::add_or_set_site_property(const std::string name,
    const double value,
    const int particle_index,
    const int site_index) {
  particles_->add_or_set_site_property(name, value,
    particle_index, site_index);
}

void Configuration::set_site_property(const std::string name,
    const double value,
    const int particle_index,
    const int site_index) {
  particles_->set_site_property(name, value, particle_index, site_index); }

void Configuration::set_site_property(const int index,
    const double value,
    const int particle_index,
    const int site_index) {
  particles_->set_site_property(index, value, particle_index, site_index); }

void Configuration::replace_properties_(const int particle_index,
                         const int site_index,
                         const Properties& prop) {
  particles_->replace_properties(particle_index, site_index, prop); }

const Particle& Configuration::particle_(const int index) {
  return particles_->particle(index); }

const std::vector<std::shared_ptr<Select> >& Configuration::ghosts() const {
  return ghosts_;
}

const std::vector<std::vector<std::shared_ptr<Table3D> > >&
  Configuration::table3d() const { return table3d_; }
const std::vector<std::vector<std::shared_ptr<Table4D> > >&
  Configuration::table4d() const { return table4d_; }
const std::vector<std::vector<std::shared_ptr<Table5D> > >&
  Configuration::table5d() const { return table5d_; }
const std::vector<std::vector<std::shared_ptr<Table6D> > >&
  Configuration::table6d() const { return table6d_; }
std::vector<std::vector<std::shared_ptr<Table3D> > > *
  Configuration::get_table3d() { return &table3d_; }
std::vector<std::vector<std::shared_ptr<Table4D> > > *
  Configuration::get_table4d() { return &table4d_; }
std::vector<std::vector<std::shared_ptr<Table5D> > > *
  Configuration::get_table5d() { return &table5d_; }
std::vector<std::vector<std::shared_ptr<Table6D> > > *
  Configuration::get_table6d() { return &table6d_; }

int Configuration::num_groups() const {
  return static_cast<int>(group_selects_.size());
}

const std::vector<std::shared_ptr<Select> >& Configuration::group_selects()
    const {
  return group_selects_;
}

const Select& Configuration::group_select(const int index) const {
  ASSERT(index < static_cast<int>(group_selects_.size()), "index:" << index
    << " >= number of groups:" << group_selects_.size());
  return *group_selects_[index];
}

const Select& Configuration::selection_of_all() const {
  return *group_selects_[0];
}

const std::string& Configuration::site_type_to_name(const int site_type) const {
  return unique_types().site_type_to_name(site_type);
}

int Configuration::site_name_to_index(const std::string& site_name) const {
  return particle_types().site_name_to_index(site_name);
}

int Configuration::site_type_name_to_index(const std::string& site_type_name) const {
  return unique_types().site_type_name_to_index(site_type_name);
}

const std::string& Configuration::site_index_to_name(const int particle_index,
    const int site_index) const {
  return particle_types().site_index_to_name(particle_index, site_index);
}

int Configuration::particle_name_to_type(const std::string& name) const {
  return unique_types_->name_to_index(name);
}

const std::string& Configuration::particle_type_to_name(const int ptype) const {
  return unique_types_->index_to_name(ptype);
}

void Configuration::site_type_names(std::vector<std::string> * names) const {
  names->resize(num_site_types());
  for (int stype = 0; stype < num_site_types(); ++stype) {
    (*names)[stype] = site_type_to_name(stype);
  }
}

}  // namespace feasst
