
#ifndef FEASST_CORE_FILE_XYZ_H_
#define FEASST_CORE_FILE_XYZ_H_

#include <fstream>
#include <sstream>
#include "core/include/utils_io.h"
#include "core/include/configuration.h"
#include "core/include/debug.h"

namespace feasst {

/// HWH Add html link to XYZ format
class FileXYZ {
 public:
  void load(const char* file_name, Configuration * config) {
    std::ifstream xyz_file(file_name);
    int num_sites;
    std::string line;
    std::getline(xyz_file, line);
    { std::stringstream iss(line);
      iss >> num_sites; }
    std::getline(xyz_file, line);
    std::vector<double> coord(3);
    { std::stringstream iss(line);
      double id;
      iss >> id >> coord[0] >> coord[1] >> coord[2]; }
    //cout << "line " << line << " is " << iss.str() << endl;
//      cout << "cord " << str(coord) << endl;
    Position position;
    position.set_vector(coord);
    DomainCuboid domain = config->domain();
    domain.set_side_length(position);
    config->set_domain(domain);
    //config->domain.set_side_length(position);
    if (config->num_sites() == 0) {
      if (config->num_particle_types() == 0) {
        Particle particle;
        particle.default_particle();
        config->add_particle_type("../forcefield/data.atom");
      }
      for (int site_index = 0; site_index < num_sites; ++site_index) {
        config->add_particle(0);
      }
    }

    /// read the coordinates into a 2d vector
    std::vector<std::vector<double> > coords;
    coords.resize(num_sites, std::vector<double>(config->dimension()));
    for (int i = 0; i < num_sites; ++i) {
      std::getline(xyz_file, line);
      std::istringstream iss(line);
      std::string tmp;
      iss >> tmp;
      for (int dim = 0; dim < config->dimension(); ++dim) {
        iss >> coords[i][dim];
      }
      DEBUG("coord " << coords[i][0]);
    }
    config->update_positions(coords);
  }

  void write(const char* file_name,
             const Configuration& config) {
    std::ofstream file(file_name);
    file << config.num_sites() << endl
         << "-1" << endl;
    // int index = 1;
    for (const Particle& part : config.particles().particles()) {
      for (const Site& site : part.sites()) {
        file << "H " << site.position().str() << endl;
      }
    }
  }
};

}  // namespace feasst

#endif  // FEASST_CORE_FILE_XYZ_H_
