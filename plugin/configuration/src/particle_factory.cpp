#include "utils/include/debug.h"
#include "utils/include/utils.h"
#include "utils/include/serialize.h"
#include "math/include/histogram.h"
#include "math/include/utils_math.h"
#include "configuration/include/file_particle.h"
#include "configuration/include/particle_factory.h"

namespace feasst {

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
  particles_.push_back(particle_copy);
  check_site_types();
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
    for (const Dihedral& dihedral : particle.dihedrals()) dihedral_type.add(dihedral.type());
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
}

void ParticleFactory::add(const std::string file_name) {
  ASSERT(unique_particles_,
    "only add particles by file for defining allowed types");
  Particle particle = FileParticle().read(file_name);

  // Assign per-site properties from the data file.
  if (unique_types_) {
    particle.remove_non_unique_types();
    FileParticle().read_properties(file_name, &particle);
  }

  add(particle);

  // Update model parameters only after the particle has been filtered.
  if (unique_types_) {
    model_params_.add(particles_.back());
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

void ParticleFactory::remove(const Group group) {
  for (int index = num() -1;
       index >= 0;
       --index) {
    Particle * part = &particles_[index];
    if (group.is_in(*part, index)) {
      group.remove_sites(part);
      //*part = group.remove_sites(*part);
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
  model_params_.serialize(ostr);
}

ParticleFactory::ParticleFactory(std::istream& istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 692, "unrecognized version: " << version);
  feasst_deserialize_fstobj(&particles_, istr);
  feasst_deserialize(&unique_particles_, istr);
  feasst_deserialize(&unique_types_, istr);
  model_params_ = ModelParams(istr);
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

}  // namespace feasst
