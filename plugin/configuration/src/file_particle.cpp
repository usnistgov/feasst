#include <string>
#include <vector>
#include <fstream>
#include "configuration/include/file_particle.h"
#include "utils/include/file.h"
#include "utils/include/debug.h"
#include "utils/include/serialize.h"

namespace feasst {

void FileParticle::read_num_and_types_(const std::string file_name) {
  std::ifstream file(file_name.c_str());
  ASSERT(file.good(), "cannot find lammps DATA file " << file_name.c_str());

  // skip all lines beginning with the character "#"
  skip_characters('#', file);


  // read number of sites
  // optionally, read number of dimentions
  std::string line, descript, descript2;
  file >> num_sites_ >> descript;
  int attempt = 0;
  const int max_attempt = 1e3;
  while (descript != "sites" && attempt < max_attempt) {
    if (descript == "dimensions") {
      num_dimensions_ = num_sites_;
    }
    file >> num_sites_ >> descript;
    ++attempt;
  }
  ASSERT(attempt != max_attempt, "unable to read number of sites in file. " <<
    "Please verify the file syntax.");

  // read next line, if it is number of bonds, then record.
  // Else, it is number of site types
  int read_num;
  file >> read_num >> descript;
  if (descript.compare("bonds") == 0) {
    num_bonds_ = read_num;

    // read next line, if it is number of angles, then record.
    // Else, it is number of site types
    file >> read_num >> descript;
    if (descript.compare("angles") == 0) {
      num_angles_ = read_num;
      // read next line, if it is number of dihedrals, then record.
      // Else, it is number of site types
      file >> read_num >> descript;
      if (descript.compare("dihedrals") == 0) {
        num_dihedrals_ = read_num;
        file >> num_site_types_ >> descript;
        ASSERT(descript.compare("site") == 0,
          "unrecognized lammps DATA format for file " << file_name);
      } else if (descript == "site") {
        num_site_types_ = read_num;
      } else {
        ASSERT(0, "unrecognized lammps DATA format for file " << file_name);
      }
    } else if (descript.compare("site") == 0) {
      num_site_types_ = read_num;
    } else {
      ASSERT(0, "unrecognized lammps DATA format for file " << file_name);
    }
  } else if (descript.compare("site") == 0) {
    num_site_types_ = read_num;
  } else {
    ASSERT(0, "unrecognized lammps DATA format for file " << file_name);
  }

  // read number of site types, if not already done so
  if (num_site_types_ == 0) {
    file >> num_site_types_ >> descript >> descript2;
    attempt = 0;
    while (descript.compare("site") != 0) {
      file >> num_site_types_ >> descript >> descript2;
      ++attempt;
    }
    ASSERT(attempt != max_attempt, "unable to read number of site types in " <<
      "file. Please verify the file syntax.");
  } else {
    file >> descript2;
  }
  ASSERT(num_site_types_ > 0, "read error");

  // read number of bond and angle types
  if (num_bonds_ != 0) {
    file >> num_bond_types_ >> descript >> descript;
  }
  if (num_angles_ != 0) {
    file >> num_angle_types_ >> descript >> descript;
  }
  if (num_dihedrals_ != 0) {
    file >> num_dihedral_types_ >> descript >> descript;
  }

  // HWH implement reading impropers
  num_impropers_ = 0;
  num_improper_types_ = 0;
}

Particle FileParticle::read(const std::string file_name) {
  feasst::Particle particle;

  std::ifstream file(file_name.c_str());
  ASSERT(file.good(), "cannot find lammps DATA file " << file_name.c_str());

  read_num_and_types_(file_name);

  // read until Sites section
  bool is_found = find("Sites", file);
  if (!is_found) FATAL("Could not find Sites in " << file_name);

  // read Sites section
  int isite, itype;
  std::vector<double> xtmp(num_dimensions_);
  std::string cm, typetmp;
  feasst::Site site;
  feasst::Position position;
  std::string line;
  for (int i = 0; i < num_sites_; ++i) {
    std::getline(file, line);
    file >> isite >> itype;
    for (int dim = 0; dim < num_dimensions_; ++dim) {
      file >> xtmp[dim];
    }

    position.set_vector(xtmp);
    site.set_position(position);
    site.set_type(itype);
    particle.add(site);
  }

  // read Bonds section
  if (num_bonds_ > 0) {
    is_found = find("Bonds", file);
    if (!is_found) FATAL("Could not find Bonds in " << file_name);
    int ibond, a1, a2;
    for (int bond_index = 0; bond_index < num_bonds_; ++bond_index) {
      file >> ibond >> itype >> a1 >> a2;
      feasst::Bond bond;
      bond.add_site_index(a1);
      bond.add_site_index(a2);
      bond.set_type(itype);
      particle.add_bond(bond);
    }
  }

  // read Angles section
  if (num_angles_ > 0) {
    is_found = find("Angles", file);
    if (!is_found) FATAL("Could not find Angles in " << file_name);
    int iangle, a1, a2, a3;
    for (int angle_index = 0; angle_index < num_angles_; ++angle_index) {
      file >> iangle >> itype >> a1 >> a2 >> a3;
      feasst::Angle angle;
      angle.add_site_index(a1);
      angle.add_site_index(a2);
      angle.add_site_index(a3);
      angle.set_type(itype);
      particle.add_angle(angle);
    }
  }

  // read Dihedrals section
  if (num_dihedrals_ > 0) {
    is_found = find("Dihedrals", file);
    if (!is_found) FATAL("Could not find Dihedrals in " << file_name);
    int idihedral, a1, a2, a3, a4;
    for (int dihedral_index = 0; dihedral_index < num_dihedrals_; ++dihedral_index) {
      file >> idihedral >> itype >> a1 >> a2 >> a3 >> a4;
      feasst::Dihedral dihedral;
      dihedral.add_site_index(a1);
      dihedral.add_site_index(a2);
      dihedral.add_site_index(a3);
      dihedral.add_site_index(a4);
      dihedral.set_type(itype);
      particle.add_dihedral(dihedral);
    }
  }

//  // read Impropers section
//  if (num_impropers_ > 0) {
//    is_found = find("Impropers", file);
//    if (!is_found) FATAL("Could not find Impropers in " << file_name);
//    int iimproper, a1, a2, a3, a4;
//    for (int improper_index = 0; improper_index < num_impropers_; ++improper_index) {
//      file >> iimproper >> itype >> a1 >> a2 >> a3 >> a4;
//      feasst::Improper improper;
//      improper.add_site_index(a1);
//      improper.add_site_index(a2);
//      improper.add_site_index(a3);
//      improper.add_site_index(a4);
//      improper.set_type(itype);
//      particle.add_improper(improper);
//    }
//  }

  return particle;
}

void FileParticle::read_properties(const std::string file_name,
                               Particle* particle) {
  std::ifstream file(file_name.c_str());
  ASSERT(file.good(), "cannot find lammps DATA file " << file_name.c_str());

  read_num_and_types_(file_name);

  bool is_found = find("Site Properties", file);
  if (!is_found) FATAL("Could not find \"Site Properties\" in " << file_name);
  read_properties_("site", num_site_types_, particle, file);
  if (num_bonds_ != 0) {
    is_found = find("Bond Properties", file);
    if (!is_found) FATAL("Could not find \"Bond Properties\" in " << file_name);
    read_properties_("bond", num_bond_types_, particle, file);
  }
  if (num_angles_ != 0) {
    is_found = find("Angle Properties", file);
    if (!is_found) FATAL("Could not find \"Angle Properties\" in " << file_name);
    read_properties_("angle", num_angle_types_, particle, file);
  }
  if (num_dihedrals_ != 0) {
    is_found = find("Dihedral Properties", file);
    if (!is_found) FATAL("Could not find \"Dihedral Properties\" in " << file_name);
    read_properties_("dihedral", num_dihedral_types_, particle, file);
  }
}

void FileParticle::read_properties_(const std::string property_type,
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
    ASSERT(properties.size() >= 3, "Missing properties for: " << property_type);
    const int type = stoi(properties[0]);
    int shift = 0;
    if (property_type == "bond") {
      shift = 1;
      particle->add_bond_model(type, properties[1]);
    } else if (property_type == "angle") {
      shift = 1;
      particle->add_angle_model(type, properties[1]);
    } else if (property_type == "dihedral") {
      shift = 1;
      particle->add_dihedral_model(type, properties[1]);
    }
    ASSERT((properties.size() - shift) % 2 == 1, "size error");
    const int num_properties = (properties.size() - 1 - shift)/2;
    DEBUG("type: " << type);
    ASSERT(type < num_types, "type: " << type << " is too large for the number "
      << "of types: " << num_types << " of property: " << property_type
      << " . Properties are listed from indices of 0 to n-1, not 1 to n. "
      << "See /path/to/feasst/forcefield/README.rst for more details.");
    for (int index = 0; index < num_properties; ++index) {
      const std::string name = properties[2*index + 1 + shift];
      const double value = stod(properties[2*index + 2 + shift]);
      if (property_type == "site") {
        Site * site = particle->get_site(type);
        DEBUG("adding " << name << " "  << value);
        site->add_property(name, value);
        if (name == "anisotropic") {
          site->set_anisotropic("true");
        }
        //particle->set_site(type, site);
      } else if (property_type == "bond") {
        DEBUG("adding bond coeff of type " << type << " name " << name << " value " << value);
        particle->add_bond_property(type, name, value);
      } else if (property_type == "angle") {
        DEBUG("adding angle coeff of type " << type << " name " << name << " value " << value);
        particle->add_angle_property(type, name, value);
      } else if (property_type == "dihedral") {
        DEBUG("adding dihedral coeff of type " << type << " name " << name << " value " << value);
        particle->add_dihedral_property(type, name, value);
//      } else if (property_type == "improper") {
//        DEBUG("adding improper coeff of type " << type << " name " << name << " value " << value);
//        particle->add_improper_property(type, name, value);
      } else {
        ERROR("unrecognized property_type: " << property_type);
      }
    }
  }
}

}  // namespace feasst
