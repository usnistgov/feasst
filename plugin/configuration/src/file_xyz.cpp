#include "utils/include/arguments.h"
#include "utils/include/io.h"
#include "utils/include/serialize.h"
#include "utils/include/debug.h"
#include "math/include/constants.h"
#include "configuration/include/particle_factory.h"
#include "configuration/include/domain.h"
#include "configuration/include/select.h"
#include "configuration/include/file_vmd.h"
#include "configuration/include/file_xyz.h"

namespace feasst {

FileXYZ::FileXYZ(argtype * args) {
  group_index_ = 0;
  if (used("group_index", *args)) {
    group_index_ = integer("group_index", args);
    ASSERT(!used("group", *args),
      "cant specify both group_index and group name");
  } else {
    if (used("group", *args)) {
      group_ = str("group", args);
    }
  }
  append_ = boolean("append", args, false);
  euler_ = boolean("euler", args, false);
}
FileXYZ::FileXYZ(argtype args) : FileXYZ(&args) {
  feasst_check_all_used(args);
}
FileXYZ::~FileXYZ() {}

bool FileXYZ::load_frame(std::ifstream& xyz_file,
                         Configuration * config) const {
  ASSERT(xyz_file, "xyz_file is empty");
  if (xyz_file.peek() == EOF) {
    return false;
  }
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
  // If third coordinate, z, is zero, then its a 2d simulation.
  if (coord[2] < NEAR_ZERO) {
    coord.pop_back();
  }
  Position position;
  position.set_vector(coord);
  DEBUG("coord " << feasst::feasst_str(coord) << " pos " << position.str());
  config->set_side_lengths(position);
  // ASSERT(config->num_sites() == num_sites, "site mismatch");
  if (config->num_sites() == 0) {
    ASSERT(config->num_particle_types() > 0, "try adding particle type " <<
      "before loading the configuration");
    ASSERT(config->num_particle_types() == 1,
      "If more than one particle type, add particles in the specific order " <<
      "before reading the XYZ file (e.g., as Configuration " <<
      "add_particles_of_type0 [n0] add_particles_of_type1 [n1] ...)");
    const int site_per_part = config->particle_types().particle(0).num_sites();
    ASSERT(num_sites % site_per_part == 0, "assumed particle type to load " <<
      "xyz file is incompatible with number of sites:" << num_sites <<
      ". site_per_part:" << site_per_part);
    for (int index = 0; index < num_sites/site_per_part; ++index) {
      config->add_particle_of_type(0);
    }
  }

  /// read the coordinates into a 2d vector
  std::vector<std::vector<double> > coords, eulers;
  coords.resize(num_sites, std::vector<double>(config->dimension()));
  eulers.resize(num_sites, std::vector<double>(config->dimension()));
  for (int i = 0; i < num_sites; ++i) {
    std::getline(xyz_file, line);
    std::istringstream iss(line);
    std::string tmp;
    iss >> tmp;
    for (int dim = 0; dim < config->dimension(); ++dim) {
      iss >> coords[i][dim];
    }
    DEBUG("coord " << coords[i][0]);
    if (euler_) {
      ASSERT(config->dimension() == 3, "Euler must have 3 dimensions.");
      for (int dim = 0; dim < config->dimension(); ++dim) {
        iss >> eulers[i][dim];
      }
      DEBUG("euler " << eulers[i][0]);
    }
  }

  // update the number of particles
  if (num_sites != config->num_sites()) {
    DEBUG("update number of particles");
    const int particle_type = 0;
    INFO("config->num_sites() " << config->num_sites());
    INFO("config->num_particle_types() " << config->num_particle_types());
    ASSERT(config->num_particle_types() == 1, "assumes 1 particle type");
    if (num_sites < config->num_sites()) {
      int spi = config->num_particles() - 1;
      Select sel;
      sel.add_particle(config->particle_type(particle_type), 0);
      while (num_sites < config->num_sites()) {
        ASSERT(spi >= 0, "error");
        sel.set_particle(0, config->selection_of_all().particle_index(spi));
        config->remove_particle(sel);
        --spi;
      }
    } else {
      while (num_sites > config->num_sites()) {
        config->add_particle_of_type(particle_type);
        ASSERT(config->num_sites() < 1e20, "infinite loop?");
      }
    }
  }

  if (euler_) {
    config->update_positions(coords, eulers);
  } else {
    config->update_positions(coords);
  }
  return true;
}

void FileXYZ::load(const std::string file_name, Configuration * config) const {
  std::ifstream xyz_file(file_name);
  ASSERT(xyz_file, "file_name: " << file_name << " is empty");
  load_frame(xyz_file, config);
}

PrinterXYZ::PrinterXYZ(std::shared_ptr<std::ofstream> file, const bool euler,
    const int num_places) : file_(file) {
  num_places_ = num_places;
  euler_ = euler;
}

void PrinterXYZ::work(const Site& site,
    const Configuration& config,
    const LoopDescriptor& data) {
  (*file_.get()) << site.type() << " ";
  (*file_.get()) << std::setprecision(num_places_);
  for (int dim = 0; dim < config.dimension(); ++dim) {
    (*file_.get()) << site.position().coord(dim) << " ";
  }
  for (int dim = config.dimension(); dim < 3; ++dim) {
    (*file_.get()) << "0 ";
  }
  if (euler_) {
    const Euler& euler = site.euler();
    (*file_.get()) << euler.phi() << " "
                   << euler.theta() << " "
                   << euler.psi();
  }
  (*file_.get()) << std::endl;
}

void FileXYZ::write(const std::string file_name,
                    const Configuration& config,
                    const int num_places) const {
  int gindex = group_index_;
  if (!group_.empty()) {
    gindex = config.group_index(group_);
  }
  auto file = std::make_shared<std::ofstream>();
  if (append_ == 0) {
    file->open(file_name);
  } else {
    file->open(file_name, std::ofstream::app);
  }
  const Domain& domain = config.domain();
  (*file.get()) << config.group_select(gindex).num_sites() << std::endl
    << "-1 ";
  (*file.get()) << std::setprecision(num_places);
  for (int dim = 0; dim < domain.dimension(); ++dim) {
    (*file.get()) << domain.side_length(dim) << " ";
  }
  (*file.get()) << domain.xy() << " "
    << domain.xz() << " "
    << domain.yz() << " "
    << std::endl;
  PrinterXYZ printer(file, euler_);
  VisitConfiguration().loop(config, &printer, gindex);
}

void FileXYZ::write_for_vmd(const std::string file_name,
    const Configuration& config) const {
  std::stringstream ss;
  ss << file_name << ".vmd";
  write(file_name, config);
  FileVMD().write(ss.str(), config, file_name);
}

void FileXYZ::serialize(std::ostream& ostr) const {
  feasst_serialize_version(2867, ostr);
  feasst_serialize(group_index_, ostr);
  feasst_serialize(group_, ostr);
  feasst_serialize(append_, ostr);
  feasst_serialize(euler_, ostr);
}

FileXYZ::FileXYZ(std::istream& istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 2867, "version mismatch: " << version);
  feasst_deserialize(&group_index_, istr);
  feasst_deserialize(&group_, istr);
  feasst_deserialize(&append_, istr);
  feasst_deserialize(&euler_, istr);
}

}  // namespace feasst
