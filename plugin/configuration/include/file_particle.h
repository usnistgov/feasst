
#ifndef FEASST_CONFIGURATION_FILE_PARTICLE_H_
#define FEASST_CONFIGURATION_FILE_PARTICLE_H_

#include <string>
#include "configuration/include/particle.h"

namespace feasst {

/**
  See <a href="../../../particle/README.html">Particle files and units</a> for more details.
 */
class FileParticle {
 public:
  /// Create a particle from the FEASST particle file format.
  Particle read(const std::string file_name);

  /// Read pair properties from file and assign properties to particle
  void read_properties(const std::string file_name, Particle* particle);

  /// Return the number of enteries in a given section, which is defined by
  /// an empty line at the beginning and end.
  /// Ignore lines beginning with a comment character.
  int read_num_in_section(const std::string& section, const std::string& file_name, const std::string comment = "#") const;

  int num_sites() const { return num_sites_; }
  int num_site_types() const { return num_site_types_; }
  int num_bonds() const { return num_bonds_; }
  int num_bond_types() const { return num_bond_types_; }
  int num_angles() const { return num_angles_; }
  int num_angle_types() const { return num_angle_types_; }
  int num_dihedrals() const { return num_dihedrals_; }
  int num_dihedral_types() const { return num_dihedral_types_; }
  int num_impropers() const { return num_impropers_; }
  int num_improper_types() const { return num_improper_types_; }

 private:
  int num_dimensions_ = 3;
  int num_sites_ = 0;
  int num_site_types_ = 0;
  int num_bonds_ = 0;
  int num_bond_types_ = 0;
  int num_angles_ = 0;
  int num_angle_types_ = 0;
  int num_dihedrals_ = 0;
  int num_dihedral_types_ = 0;
  int num_impropers_ = 0;
  int num_improper_types_ = 0;

  void read_num_and_types_(const std::string file_name);

  void read_properties_(const std::string property_type,
                        const int num_types,
                        Particle * particle,
                        std::ifstream & file) const;

  // return the type as an index given type name and list of types
  int find_type_index_(const std::string& type,
    std::vector<std::string> * types) const;

  // return the index given the name
  int name_to_index(const std::string& name, const std::vector<std::string>& names) const;
};

}  // namespace feasst

#endif  // FEASST_CONFIGURATION_FILE_PARTICLE_H_
