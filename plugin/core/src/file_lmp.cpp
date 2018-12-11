#include <string>
#include <vector>
#include <fstream>
#include "core/include/file_lmp.h"
#include "core/include/file.h"
#include "core/include/debug.h"
#include "core/include/utils_io.h"

namespace feasst {

void FileLMP::read_num_and_types_(const std::string file_name) {
  // open LAMMPS data file
  std::ifstream file(file_name.c_str());
  ASSERT(file.good(), "cannot find lammps DATA file " << file_name.c_str());

  // skip all lines beginning with the character "#"
  skip_characters('#', file);

  // read number of atoms
  std::string line, descript, descript2;
  file >> num_atoms_ >> descript;
  while (descript.compare("atoms") != 0) {
    file >> num_atoms_ >> descript;
  }

  // read next line, if it is number of bonds, then record.
  // Else, it is number of atom types
  int read_num;
  file >> read_num >> descript;
  if (descript.compare("bonds") == 0) {
    num_bonds_ = read_num;

    // read next line, if it is number of angles, then record.
    // Else, it is number of atom types
    file >> read_num >> descript;
    if (descript.compare("angles") == 0) {
      num_angles_ = read_num;
      file >> num_atom_types_ >> descript;
      ASSERT(descript.compare("atom") == 0,
        "unrecognized lammps DATA format for file " << file_name);
    } else if (descript.compare("atom") == 0) {
      num_atom_types_ = read_num;
    } else {
      ASSERT(0, "unrecognized lammps DATA format for file " << file_name);
    }
  } else if (descript.compare("atom") == 0) {
    num_atom_types_ = read_num;
  } else {
    ASSERT(0, "unrecognized lammps DATA format for file " << file_name);
  }

  // read number of atom types, if not already done so
  if (num_atom_types_ == 0) {
    file >> num_atom_types_ >> descript >> descript2;
    while (descript.compare("atom") != 0) {
      file >> num_atom_types_ >> descript >> descript2;
    }
  } else {
    file >> descript2;
  }
  ASSERT(num_atom_types_ > 0, "read error");

  // read number of bond and angle types
  if (num_bonds_ != 0) {
    file >> num_bond_types_ >> descript >> descript;
  }
  if (num_angles_ != 0) {
    file >> num_angle_types_ >> descript >> descript;
  }
}

Particle FileLMP::read(const std::string file_name) {
  feasst::Particle particle;

  // open LAMMPS data file
  std::ifstream file(file_name.c_str());
  ASSERT(file.good(), "cannot find lammps DATA file " << file_name.c_str());

  read_num_and_types_(file_name);

  // read until Atoms section
  find_or_fail("Atoms", file);

  // read Atoms section
  int iatom, imol, itype;
  const int dimension = 3;
  std::vector<double> xtmp(dimension);
  std::string cm, typetmp;
  feasst::Site site;
  feasst::Position position;
  std::string line;
  for (int i = 0; i < num_atoms_; ++i) {
    std::getline(file, line);
    file >> iatom >> imol >> itype;
    for (int dim = 0; dim < dimension; ++dim) {
      file >> xtmp[dim];
    }

    position.set_vector(xtmp);
    if (i == 0) {
      particle.set_position(position);
    }
    site.set_position(position);
    site.set_type(itype - 1);
    particle.add(site);
  }

  // read Bonds section
  if (num_bonds_ > 0) {
    find_or_fail("Bonds", file);
    int ibond, a1, a2;
    for (int bond_index = 0; bond_index < num_bonds_; ++bond_index) {
      file >> ibond >> itype >> a1 >> a2;
      feasst::Bond bond;
      bond.add_site_index(a1 - 1);
      bond.add_site_index(a2 - 1);
      bond.set_type(itype - 1);
      particle.add_bond(bond);
    }
  }

  return particle;
}

void FileLMP::read_properties(const std::string file_name,
                               Particle* particle) {
  // open LAMMPS data file
  std::ifstream file(file_name.c_str());
  ASSERT(file.good(), "cannot find lammps DATA file " << file_name.c_str());

  read_num_and_types_(file_name);

  find_or_fail("Atom Properties", file);
  read_properties_("site", num_atom_types_, particle, file);
  if (num_bonds_ != 0) {
    find_or_fail("Bond Properties", file);
    read_properties_("bond", num_bond_types_, particle, file);
  }
}

void FileLMP::read_properties_(const std::string property_type,
                               const int num_types,
                               Particle * particle,
                               std::ifstream & file) const {
  std::string line;
  std::getline(file, line);
  DEBUG("read properties: " << line);
  for (int i = 0; i < num_types; ++i) {
    std::getline(file, line);
    DEBUG("read properties i " << i << ": " << line);
    std::vector<std::string> properties = split(line);
    ASSERT(properties.size() >= 3, "size error");
    ASSERT(properties.size() % 2 == 1, "size error");
    const int num_properties = (properties.size() - 1)/2;
    const int type = stoi(properties[0]) - 1;
    for (int index = 0; index < num_properties; ++index) {
      const std::string name = properties[2*index + 1];
      const double value = stod(properties[2*index + 2]);
      if (property_type == "site") {
        Site site = particle->site(type);
        site.add_property(name, value);
        particle->set_site(type, site);
      } else if (property_type == "bond") {
        DEBUG("adding bond coeff of type " << type << " name " << name << " value " << value);
        Bond bond = particle->bond(type);
        bond.add_property(name, value);
        particle->set_bond(type, bond);
      } else {
        ERROR("unrecognized property_type: " << property_type);
      }
    }
  }
}

}  // namespace feasst
