#include <string>
#include <iostream>
#include "utils/include/debug.h"
#include "utils/include/utils.h"
#include "utils/include/serialize.h"
#include "math/include/histogram.h"
#include "math/include/utils_math.h"
#include "configuration/include/model_params.h"
#include "configuration/include/group.h"
#include "configuration/include/file_particle.h"
#include "configuration/include/particle_factory.h"

namespace feasst {

ParticleFactory::ParticleFactory() {
  model_params_ = std::make_shared<ModelParams>();
}

void ParticleFactory::add(const Particle& particle) {
  // compute types and also check none are skipped.
  const int num_site_type = num_site_types();
  const int num_particle_type = num_particle_types();
  Particle particle_copy = particle;
  if (unique_types_) {
    particle_copy.remove_non_unique_types();
  }
  if (unique_particles_) {
    particle_copy.increment_site_types(num_site_type);
    particle_copy.set_type(num_particle_type);
  }
  if (!unique_types_ && !unique_particles_) {
    particle_copy.clear_names();
  }
  particles_.push_back(particle_copy);
  if (unique_types_ || unique_particles_) {
    check_site_types();
  }
}

void ParticleFactory::check(const int dimension) const {
  if (particles_.size() != 0) {
    int size = dimension;
    if (size == -1) {
      size = particles_.front().site(0).position().size();
    }
    for (Particle particle : particles_) {
      ASSERT(particle.site(0).position().size() == size, "size error");
      particle.check();
    }
  }
}

void ParticleFactory::check_types(int * num_site_types,
    int * num_particle_types, int * num_bond_types, int * num_angle_types,
    int * num_dihedral_types) const {
  if (particles_.size() == 0) {
    *num_site_types = 0;
    *num_particle_types = 0;
    *num_bond_types = 0;
    *num_angle_types = 0;
    *num_dihedral_types = 0;
    return;
  }
  if (!unique_types_ && !unique_particles_) {
    return;
  }

  Histogram site_type;
  site_type.set_width_center(1, 0);
  Histogram particle_type(site_type), bond_type(site_type),
    angle_type(site_type), dihedral_type(site_type);
  for (const Particle& particle : particles_) {
    TRACE("particle type " << particle.type());
    particle_type.add(particle.type());
    for (const Site& site : particle.sites()) site_type.add(site.type());
    for (const Bond& bond : particle.bonds()) bond_type.add(bond.type());
    for (const Angle& angle : particle.angles()) angle_type.add(angle.type());
    for (const Dihedral& dihedral : particle.dihedrals()) {
      dihedral_type.add(dihedral.type());
    }
  }
  if (unique_particles_) {
    for (const double& value : site_type.histogram()) {
      ASSERT(value > 0, "particle skipped a site type");
    }
    for (const double& value : particle_type.histogram()) {
      ASSERT(value > 0, "skipped a particle type");
    }
    if (sum(bond_type.histogram()) > 0) {
      for (const double& value : bond_type.histogram()) {
        ASSERT(value > 0, "particle skipped a bond type");
      }
    }
    if (sum(angle_type.histogram()) > 0) {
      for (const double& value : angle_type.histogram()) {
        ASSERT(value > 0, "particle skipped an angle type");
      }
    }
    if (sum(dihedral_type.histogram()) > 0) {
      for (const double& value : dihedral_type.histogram()) {
        ASSERT(value > 0, "particle skipped an dihedral type");
      }
    }
  }
  *num_site_types = round(site_type.center_of_last_bin()) + 1;
  *num_particle_types = round(particle_type.center_of_last_bin()) + 1;
  *num_bond_types = round(bond_type.center_of_last_bin()) + 1;
  *num_angle_types = round(angle_type.center_of_last_bin()) + 1;
  *num_dihedral_types = round(dihedral_type.center_of_last_bin()) + 1;

//  // Check that all site names are unique
//  if (unique_types_) {
//    std::vector<std::string> names;
//    for (const Particle& particle : particles_) {
//      for (const Site& site : particle.sites()) {
//        ASSERT(!find_in_list(site.name(), names), "A site_type with name: " <<
//          site.name() << " is equivalent to a previous site_type name.");
//        names.push_back(site.name());
//      }
//    }
//  }
}

void ParticleFactory::rename_nonunique_(const std::vector<std::string>& names,
    const std::string& original_name, const std::string& prop, const std::string& file_name,
    const int max_warn, int * num_warn, std::string * name) {
  DEBUG("name " << *name << " names " << feasst_str(names));
  int num_attempts = 0;
  while (find_in_list(*name, names)) {
    DEBUG("renaming " << *name << " to: " << *name << "-");
    if (num_attempts > 1e2) {
      *name = *name + "-"; // faster than stoi + to_string
    } else {
      try {
        int index = std::stoi(*name);
        *name = std::to_string(index + 1);
      } catch (...) {
        *name = *name + "-";
      }
    }
    ASSERT(num_attempts < 1e6, "error");
    ++num_attempts;
    if (*num_warn < max_warn && !find_in_list(*name, names)) {
      std::string type = "name";
      if (unique_types_) type = "type";
      std::cout << "# Renamed " << prop << " " << type << ":" <<
        original_name << "->" << *name << " for particle_type:" << num() <<
        " in:" << file_name << std::endl;
     *num_warn += 1;
    }
    if (*num_warn == max_warn) {
      std::cout << "# ... and so on (suppressed following rename warnings)" << std::endl;
      *num_warn += 1;
    }
  }
}

void ParticleFactory::add(const std::string file_name) {
  ASSERT(unique_particles_,
    "only add particles by file for defining allowed types");
  Particle particle = FileParticle().read(file_name);

  DEBUG("Assign per-site properties from the data file:" << file_name);
  if (unique_types_) {
    particle.remove_non_unique_types();
    FileParticle().read_properties(file_name, &particle);
  }

  add(particle);

  if (unique_types_ || unique_particles_) {
    DEBUG("Ensure all site, bond, angle and dihedral names are unique");
    int num_warn = 0, max_warn = 2e0;
    std::vector<std::string> snames, bnames, anames, dnames;
    for (Particle& particle : particles_) {
      DEBUG("particle.type() " << particle.type());
      for (int isite = 0; isite < particle.num_sites(); ++isite) {
        Site * site = particle.get_site(isite);
        std::string name = site->name(), original_name = name;
        rename_nonunique_(snames, original_name, "Site", file_name, max_warn, &num_warn, &name);
        site->set_name(name);
        snames.push_back(site->name());
      }
      for (int ibond = 0; ibond < particle.num_bonds(); ++ibond) {
        Bond * bond = particle.get_bond(ibond);
        std::string name = bond->name(), original_name = name;
        rename_nonunique_(bnames, original_name, "Bond", file_name, max_warn, &num_warn, &name);
        bond->set_name(name);
        bnames.push_back(bond->name());
      }
      for (int iangle = 0; iangle < particle.num_angles(); ++iangle) {
        Angle * angle = particle.get_angle(iangle);
        std::string name = angle->name(), original_name = name;
        rename_nonunique_(anames, original_name, "Angle", file_name, max_warn, &num_warn, &name);
        angle->set_name(name);
        anames.push_back(angle->name());
      }
      for (int idihedral = 0; idihedral < particle.num_dihedrals(); ++idihedral) {
        Dihedral * dihedral = particle.get_dihedral(idihedral);
        std::string name = dihedral->name(), original_name = name;
        rename_nonunique_(dnames, original_name, "Dihedral", file_name, max_warn, &num_warn, &name);
        dihedral->set_name(name);
        dnames.push_back(dihedral->name());
      }
    }
  }

  DEBUG("Update model parameters only after the particle has been filtered.");
  if (unique_types_) {
    model_params_->add(particles_.back());
  }
}

ParticleFactory& ParticleFactory::unique_particles() {
  unique_particles_ = true;
  check_types();
  ASSERT(num_site_types() == 0, "set before adding particles");
  return *this;
}

ParticleFactory& ParticleFactory::unique_types() {
  unique_particles();  // one of each requires unique site types.
  unique_types_ = true;
  ASSERT(num_site_types() == 0, "set before adding particles");
  return *this;
}

int ParticleFactory::num_sites() const {
  int num = 0;
  for (const Particle& particle : particles_) {
    num += particle.num_sites();
  }
  return num;
}

int ParticleFactory::num_bonds() const {
  int num = 0;
  for (const Particle& particle : particles_) {
    num += particle.num_bonds();
  }
  return num;
}

int ParticleFactory::num_angles() const {
  int num = 0;
  for (const Particle& particle : particles_) {
    num += particle.num_angles();
  }
  return num;
}

int ParticleFactory::num_dihedrals() const {
  int num = 0;
  for (const Particle& particle : particles_) {
    num += particle.num_dihedrals();
  }
  return num;
}

void ParticleFactory::check_types() const {
  int num_site_types, num_particle_types, num_bond_types, num_angle_types,
      num_dihedral_types;
  check_types(&num_site_types, &num_particle_types,
              &num_bond_types, &num_angle_types, &num_dihedral_types);
}

int ParticleFactory::check_site_types() const {
  int num_site_types, num_particle_types, num_bond_types, num_angle_types,
      num_dihedral_types;
  check_types(&num_site_types, &num_particle_types,
              &num_bond_types, &num_angle_types, &num_dihedral_types);
  return num_site_types;
}

int ParticleFactory::check_particle_types() const {
  int num_site_types, num_particle_types, num_bond_types, num_angle_types,
      num_dihedral_types;
  check_types(&num_site_types, &num_particle_types,
              &num_bond_types, &num_angle_types, &num_dihedral_types);
  return num_particle_types;
}

int ParticleFactory::check_bond_types() const {
  int num_site_types, num_particle_types, num_bond_types, num_angle_types,
      num_dihedral_types;
  check_types(&num_site_types, &num_particle_types,
              &num_bond_types, &num_angle_types, &num_dihedral_types);
  return num_bond_types;
}

int ParticleFactory::check_angle_types() const {
  int num_site_types, num_particle_types, num_bond_types, num_angle_types,
      num_dihedral_types;
  check_types(&num_site_types, &num_particle_types,
              &num_bond_types, &num_angle_types, &num_dihedral_types);
  return num_angle_types;
}

int ParticleFactory::check_dihedral_types() const {
  int num_site_types, num_particle_types, num_bond_types, num_angle_types,
      num_dihedral_types;
  check_types(&num_site_types, &num_particle_types,
              &num_bond_types, &num_angle_types, &num_dihedral_types);
  return num_dihedral_types;
}

void ParticleFactory::remove(const Group& group) {
  for (int index = num() -1;
       index >= 0;
       --index) {
    Particle * part = &particles_[index];
    if (group.is_in(*part, index)) {
      group.remove_sites(part);
      // *part = group.remove_sites(*part);
    } else {
      particles_.erase(particles_.begin() + index);
    }
  }
}

void ParticleFactory::replace_position(const int particle_index,
                                 const Particle& replacement) {
  particles_[particle_index].replace_position(replacement);
}

void ParticleFactory::replace_position(const int particle_index,
                                 const int site_index,
                                 const Position& replacement) {
  particles_[particle_index].replace_position(site_index, replacement);
}

void ParticleFactory::remove(const int particle_index) {
  ASSERT(particle_index < num(), "size error");
  particles_.erase(particles_.begin() + particle_index);
}

void ParticleFactory::serialize(std::ostream& ostr) const {
  feasst_serialize_version(692, ostr);
  feasst_serialize_fstobj(particles_, ostr);
  feasst_serialize(unique_particles_, ostr);
  feasst_serialize(unique_types_, ostr);
  feasst_serialize(model_params_, ostr);
}

ParticleFactory::ParticleFactory(std::istream& istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 692, "unrecognized version: " << version);
  feasst_deserialize_fstobj(&particles_, istr);
  feasst_deserialize(&unique_particles_, istr);
  feasst_deserialize(&unique_types_, istr);
// HWH for unknown reasons, this function template does not work.
//  feasst_deserialize(model_params_, istr);
  { int existing;
    istr >> existing;
    if (existing != 0) {
      model_params_ = std::make_shared<ModelParams>(istr);
    }
  }
}

void ParticleFactory::scale_particle_positions(const int dimen,
    const double factor) {
  Position displacement;
  for (Particle& particle : particles_) {
    if (dimen == -1) {
      displacement = particle.site(0).position();
      displacement.multiply(factor - 1.);
    } else {
      displacement.set_to_origin(particle.site(0).position().dimension());
      displacement.set_coord(dimen,
        (factor - 1)*particle.site(0).position().coord(dimen));
    }
    particle.displace(displacement);
  }
}

void ParticleFactory::set_site_physical(const int particle, const int site,
    const bool phys) {
  ASSERT(particle < num(), "particle:" << particle << " >= number of particles:"
    << num());
  particles_[particle].set_site_physical(site, phys);
}

void ParticleFactory::replace_properties(const int particle_index,
    const int site_index,
    const Properties& replacement) {
  ASSERT(particle_index < num(), "particle_index:" << particle_index <<
    " is >= number of particles:" << num());
  particles_[particle_index].replace_properties(site_index, replacement);
}

const ModelParams& ParticleFactory::model_params() const {
  return *model_params_;
}

void ParticleFactory::add(const std::shared_ptr<ModelParam> param) {
  model_params_->add(param);
}

void ParticleFactory::set_model_param(const std::string name,
                     const int site_type,
                     const double value) {
  model_params_->set(name, site_type, value);
}

void ParticleFactory::set_model_param(const std::string name,
                     const int site_type1,
                     const int site_type2,
                     const double value) {
  model_params_->set(name, site_type1, site_type2, value);
}

void ParticleFactory::set_model_param(const std::string name,
                                      const std::string filename) {
  model_params_->set(name, filename);
}

void ParticleFactory::add_model_param(const std::string name,
                     const double value) {
  model_params_->add_property(name, value);
}

void ParticleFactory::add_or_set_model_param(const std::string name,
                            const double value) {
  model_params_->add_or_set_property(name, value);
}

void ParticleFactory::set_cutoff_min_to_sigma() {
  model_params_->set_cutoff_min_to_sigma();
}

void ParticleFactory::set_physical_constants(
    std::shared_ptr<PhysicalConstants> constants) {
  model_params_->set_physical_constants(constants);
}

int ParticleFactory::site_name_to_index(const std::string& site_name,
    int * particle_index) const {
  ASSERT(unique_particles_, "only logical to use for unique_particles");
  for (int pindex = 0; pindex < num(); ++pindex) {
    const Particle& part = particles_[pindex];
    for (int site_index = 0; site_index < part.num_sites(); ++site_index) {
      const Site& site = part.site(site_index);
      if (site.name() == site_name) {
        if (particle_index) *particle_index = pindex;
        return site_index;
      }
    }
  }
  FATAL("Could not find a site_name: " << site_name);
}

int ParticleFactory::site_type_to_particle_type(const int site_type) const {
  ASSERT(unique_types_, "only logical to use for unique_types");
  int particle_type = 0;
  int prev = 0;
  for (int site = 0; site <= site_type; ++site) {
    const int num_sites = particle(particle_type).num_sites();
    if (site >= num_sites + prev) {
      prev += num_sites;
      ++particle_type;
    }
  }
  return particle_type;
}

const Site& ParticleFactory::unique_type(const int ptype,
    const int stype) const {
  ASSERT(unique_types_, "only logical to use for unique_types");
  int prev_site_types = 0;
  for (int pt = 0; pt < ptype; ++pt) {
    prev_site_types += particle(pt).num_sites();
  }
  const int index = stype - prev_site_types;
  ASSERT(index >= 0, "error");
  return particle(ptype).site(index);
}

const std::string& ParticleFactory::site_type_to_name(const int site_type) const {
  const int ptype = site_type_to_particle_type(site_type);
  return unique_type(ptype, site_type).name();
}

int ParticleFactory::site_type_name_to_index(const std::string& site_type_name) const {
  ASSERT(unique_types_, "only logical to use for unique_types");
  int index = 0;
  for (const Particle& part : particles_) {
    for (const Site& site : part.sites()) {
      if (site.name() == site_type_name) {
        return index;
      } else {
        ++index;
      }
    }
  }
  FATAL("site_type_name:" << site_type_name << " not found.");
}

}  // namespace feasst
