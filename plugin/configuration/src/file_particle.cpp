#include <string>
#include <vector>
#include <fstream>
#include "utils/include/utils.h" // find_in_list
#include "utils/include/arguments.h"
#include "utils/include/file.h"
#include "utils/include/debug.h"
#include "utils/include/io.h"
#include "configuration/include/file_particle.h"

namespace feasst {

int FileParticle::read_section(const std::string& section,
    const std::string& file_name,
    const std::string comment,
    std::vector<std::string> * names) const {
  std::ifstream file(file_name.c_str());
  ASSERT(file.good(), "cannot find file " << file_name.c_str());
  if (!find(section, file)) {
    return 0;
  }
  std::string line;
  std::getline(file, line);
  ASSERT(line.empty(), "The \"" << section << "\" section in file:" <<
    file_name << " must begin with an empty line. Instead, found: " << line);
  std::getline(file, line);
  ASSERT(!line.empty(), "The \"" << section << "\" section in file:" <<
    file_name << " must begin with only one empty line.");
  DEBUG("line: " << line << " eof? " << file.eof());
  if (names) names->clear();
  int num = 0;
  bool last_line = false;
  while (!line.empty() && !last_line) {
    if (line.front() != comment.front()) {
      if (names && !line.empty()) {
        DEBUG("line: " << line);
        names->push_back(split(line, ' ')[0]);
      }
    }
    if (file.eof()) last_line = true;
    std::getline(file, line);
    DEBUG("line: " << line << " eof? " << file.eof());
    if (line.front() != comment.front()) {
      ++num;
    }
    DEBUG("num: " << num);
    ASSERT(num < 1e8, "The \"" << section << "\" section in file:" <<
      file_name << " must end with an empty line.");
  }
  DEBUG("num: " << num << " in Section: " << section);
  return num;
}

void read_name_in_section_(const std::string& section,
                           std::vector<std::string> * names) {
}

void FileParticle::read_num_and_types_(const std::string file_name) {
  std::ifstream file(file_name.c_str());
  ASSERT(file.good(), "cannot find file " << file_name.c_str());
  skip_characters('#', file);
  if (find("2 dimensions", file)) {
    num_dimensions_ = 2;
  }
  num_sites_ = read_section("Sites", file_name);
  num_site_types_ = read_section("Site Properties", file_name, "#", &stypes_);
  DEBUG("read site types: " << feasst_str(stypes_));
  num_bond_types_ = read_section("Bond Properties", file_name, "#", &btypes_);
  DEBUG("read bond types: " << feasst_str(btypes_));
  num_angle_types_ = read_section("Angle Properties", file_name, "#", &atypes_);
  DEBUG("read angle types: " << feasst_str(atypes_));
  num_dihedral_types_ = read_section("Dihedral Properties", file_name, "#", &dtypes_);
  DEBUG("read dihdl types: " << feasst_str(dtypes_));
  num_bonds_ = read_section("Bonds", file_name);
  num_angles_ = read_section("Angles", file_name);
  num_dihedrals_ = read_section("Dihedrals", file_name);
  num_improper_types_ = read_section("Improper Properties", file_name);
  ASSERT(num_improper_types_ == 0, "Impropers not implemented");
  num_impropers_ = read_section("Impropers", file_name);
  ASSERT(num_impropers_ == 0, "Impropers not implemented");
}

int FileParticle::find_type_index_(const std::string& type,
    const std::vector<std::string>& types) const {
  int type_index;
  //if (!find_in_list(type, *types, &type_index)) {
  //  FATAL("type:" << type << " is not found in types:" << feasst_str(*types));
  //  types->push_back(type);
  //}
  const bool found = find_in_list(type, types, &type_index);
  ASSERT(found, "type:" << type << " is not found in types:" << feasst_str(types));
  DEBUG("type_index:" << type_index);
  return type_index;
}

int FileParticle::name_to_index(const std::string& name, const std::vector<std::string>& names) const {
  int index;
  const bool found = find_in_list(name, names, &index);
  ASSERT(found, "Could not find: " << name << " in: " << feasst_str(names));
  return index;
}

Particle FileParticle::read(const std::string file_name) {
  feasst::Particle particle;

  std::ifstream file(file_name.c_str());
  ASSERT(file.good(), "cannot find file " << file_name.c_str());

  read_num_and_types_(file_name);

  // read until Sites section
  bool is_found = find("Sites", file);
  if (!is_found) FATAL("Could not find Sites in " << file_name);

  // read Sites section
  std::string sname, stype, btype, atype, dtype, al1, al2, al3, al4;
  int a1, a2, a3, a4;
  std::vector<std::string> snames(num_sites_), stypes, btypes, atypes, dtypes;
  std::vector<double> xtmp(num_dimensions_);
  std::string cm, typetmp;
  feasst::Site site;
  feasst::Position position;
  std::string line;
  for (int site_index = 0; site_index < num_sites_; ++site_index) {
    std::getline(file, line);
    file >> sname >> stype;
    for (int dim = 0; dim < num_dimensions_; ++dim) {
      file >> xtmp[dim];
    }
    position.set_vector(xtmp);
    site.set_position(position);
    const int type = find_type_index_(stype, stypes_);
    site.set_type(type);
    site.set_name(sname);
    DEBUG("stype: " << stype << " sname: " << sname << " type: " << type);
    snames[site_index] = sname;
    particle.add(site);
  }
  DEBUG("snames:" << feasst_str(snames));
  DEBUG("stypes:" << feasst_str(stypes_));

  // read Bonds section
  if (num_bonds_ > 0) {
    is_found = find("Bonds", file);
    if (!is_found) FATAL("Could not find Bonds in " << file_name);
    std::string bname;
    for (int bond_index = 0; bond_index < num_bonds_; ++bond_index) {
      file >> bname >> btype >> al1 >> al2;
      feasst::Bond bond;
      a1 = name_to_index(al1, snames);
      a2 = name_to_index(al2, snames);
      bond.add_site_index(a1);
      bond.add_site_index(a2);
      bond.set_type(find_type_index_(btype, btypes_));
      bond.set_name(bname);
      particle.add_bond(bond);
    }
  }

  // read Angles section
  if (num_angles_ > 0) {
    is_found = find("Angles", file);
    if (!is_found) FATAL("Could not find Angles in " << file_name);
    std::string aname;
    for (int angle_index = 0; angle_index < num_angles_; ++angle_index) {
      file >> aname >> atype >> al1 >> al2 >> al3;
      feasst::Angle angle;
      a1 = name_to_index(al1, snames);
      a2 = name_to_index(al2, snames);
      a3 = name_to_index(al3, snames);
      angle.add_site_index(a1);
      angle.add_site_index(a2);
      angle.add_site_index(a3);
      angle.set_type(find_type_index_(atype, atypes_));
      angle.set_name(aname);
      particle.add_angle(angle);
    }
  }

  // read Dihedrals section
  if (num_dihedrals_ > 0) {
    is_found = find("Dihedrals", file);
    if (!is_found) FATAL("Could not find Dihedrals in " << file_name);
    std::string dname;
    for (int dihedral_index = 0; dihedral_index < num_dihedrals_;
         ++dihedral_index) {
      file >> dname >> dtype >> al1 >> al2 >> al3 >> al4;
      feasst::Dihedral dihedral;
      a1 = name_to_index(al1, snames);
      a2 = name_to_index(al2, snames);
      a3 = name_to_index(al3, snames);
      a4 = name_to_index(al4, snames);
      dihedral.add_site_index(a1);
      dihedral.add_site_index(a2);
      dihedral.add_site_index(a3);
      dihedral.add_site_index(a4);
      dihedral.set_type(find_type_index_(dtype, dtypes_));
      dihedral.set_name(dname);
      particle.add_dihedral(dihedral);
    }
  }

// read Impropers section
  return particle;
}

void FileParticle::read_properties(const std::string file_name,
                               Particle* particle) {
  DEBUG("file_name: " << file_name);
  std::ifstream file(file_name.c_str());
  ASSERT(file.good(), "cannot find file " << file_name.c_str());

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
    if (!is_found) {
      FATAL("Could not find \"Angle Properties\" in " << file_name);
    }
    read_properties_("angle", num_angle_types_, particle, file);
  }
  if (num_dihedrals_ != 0) {
    is_found = find("Dihedral Properties", file);
    if (!is_found) {
      FATAL("Could not find \"Dihedral Properties\" in " << file_name);
    }
    read_properties_("dihedral", num_dihedral_types_, particle, file);
  }
}

void FileParticle::read_properties_(const std::string property_type,
    const int num_types, Particle * particle, std::ifstream & file) const {
  std::string line;
  std::getline(file, line);
  DEBUG("property_type " << property_type);
  DEBUG("read properties: " << line);
  for (int i = 0; i < num_types; ++i) {
    int shift = 0;
    if (property_type == "bond" || property_type == "angle" ||
        property_type == "dihedral") {
      shift = 1;
    }
    std::getline(file, line);
    DEBUG("read properties i " << i << ": " << line);
    std::vector<std::string> properties;
    if (is_found_in(line, "=")) {
      DEBUG("Parse new delimitor with = sign for assigning values");
      std::stringstream ss(line);
      std::string type;
      ss >> type;
      properties.push_back(type);
      for (int ishift = 0; ishift < shift; ++ishift) {
        ss >> type;
        properties.push_back(type);
      }
      std::string remaining;
      std::getline(ss, remaining);
      std::vector<std::string> pairs = split(remaining, ' ');
      for (const std::string& pair : pairs) {
        std::vector<std::string> vals = split(pair, '=');
        DEBUG("vals:" << feasst_str(vals));
        properties.push_back(vals[0]);
        properties.push_back(vals[1]);
      }
    } else {
      properties = split(line, ' ');
    }
    ASSERT(properties.size() >= 3, "Missing properties for: " << property_type);
    std::vector<std::string> types;
    if (property_type == "site") {
      types = stypes_;
    } else if (property_type == "bond") {
      types = btypes_;
    } else if (property_type == "angle") {
      types = atypes_;
    } else if (property_type == "dihedral") {
      types = dtypes_;
    } else {
      FATAL("Unrecognized property_type:" << property_type);
    }
    DEBUG("properties[0]:" << properties[0]);
    const int type = find_type_index_(properties[0], types);
    DEBUG("property_type: " << property_type);
    DEBUG("type: " << type);
    DEBUG("types: " << feasst_str(types));
    ASSERT(type < num_types, "type: " << type << " is too large for the number "
      << "of types: " << num_types << " of property: " << property_type
      << ". See /path/to/feasst/particle/README.rst for more details.");
    if (property_type == "bond") {
      particle->add_bond_model(type, properties[1]);
    } else if (property_type == "angle") {
      particle->add_angle_model(type, properties[1]);
    } else if (property_type == "dihedral") {
      particle->add_dihedral_model(type, properties[1]);
    }
    ASSERT((properties.size() - shift) % 2 == 1, "size error");
    const int num_properties = (properties.size() - 1 - shift)/2;
    for (int index = 0; index < num_properties; ++index) {
      const std::string name = properties[2*index + 1 + shift];
      const double value = stod(properties[2*index + 2 + shift]);
      if (property_type == "site") {
        Site * site = particle->get_site(type);
        if (index == 0) {
          site->set_name(properties[0]);
        }
        DEBUG("adding " << name << " "  << value);
        site->add_property(name, value);
        if (name == "anisotropic") {
          site->set_anisotropic("true");
        }
        // particle->set_site(type, site);
      } else if (property_type == "bond") {
        DEBUG("adding bond coeff of type " << type << " name " << name <<
              " value " << value);
        if (index == 0) {
          particle->get_bond(type)->set_name(properties[0]);
        }
        particle->add_bond_property(type, name, value);
      } else if (property_type == "angle") {
        DEBUG("adding angle coeff of type " << type << " name " << name <<
              " value " << value);
        if (index == 0) {
          particle->get_angle(type)->set_name(properties[0]);
        }
        particle->add_angle_property(type, name, value);
      } else if (property_type == "dihedral") {
        DEBUG("adding dihedral coeff of type " << type << " name " << name <<
              " value " << value);
        if (index == 0) {
          particle->get_dihedral(type)->set_name(properties[0]);
        }
        particle->add_dihedral_property(type, name, value);
//      } else if (property_type == "improper") {
//        particle->add_improper_property(type, name, value);
      } else {
        ERROR("unrecognized property_type: " << property_type);
      }
    }
  }
}

}  // namespace feasst
