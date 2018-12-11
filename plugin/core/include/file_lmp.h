
#ifndef FEASST_CORE_FILE_LMP_H_
#define FEASST_CORE_FILE_LMP_H_

#include "core/include/particle.h"

namespace feasst {

/// HWH link LMP format description
/// HWH: some major modifications: the charge is read in Pair Coeffs and no wrapping tags
/// HWH: For now its hard coded for 3 pair properties
class FileLMP {
 public:
  /// Create a particle from LAMMPS file format.
  Particle read(const std::string file_name);

  /// Read pair properties from file and assign properties to particle
  void read_properties(const std::string file_name, Particle* particle);

  int num_atoms() const { return num_atoms_; }
  int num_atom_types() const { return num_atom_types_; }
  int num_bonds() const { return num_bonds_; }
  int num_bond_types() const { return num_bond_types_; }
  int num_angles() const { return num_angles_; }
  int num_angle_types() const { return num_angle_types_; }

 private:
  int num_atoms_ = 0;
  int num_atom_types_ = 0;
  int num_bonds_ = 0;
  int num_bond_types_ = 0;
  int num_angles_ = 0;
  int num_angle_types_ = 0;

  void read_num_and_types_(const std::string file_name);

  void read_properties_(const std::string property_type,
                        const int num_types,
                        Particle * particle,
                        std::ifstream & file) const ;
};

}  // namespace feasst

#endif  // FEASST_CORE_FILE_LMP_H_
