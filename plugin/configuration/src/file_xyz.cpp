
#include <fstream>
#include <sstream>
#include "configuration/include/file_xyz.h"
#include "configuration/include/domain.h"
#include "utils/include/serialize.h"
#include "utils/include/debug.h"

namespace feasst {

void FileVMD::serialize(std::ostream& ostr) const {
  feasst_serialize_version(1, ostr);
}

FileVMD::FileVMD(std::istream& istr) {
  feasst_deserialize_version(istr);
}

void FileVMD::get_params(const Configuration& config,
    const int site_type,
    double * radius,
    double * distance,
    int * center_index) const {
  *radius = 0.5*config.model_params().select("sigma").value(site_type);
  *distance = 0.;
  *center_index = -1;
}

void FileVMD::write(const std::string file_name,
    const Configuration& config,
    const std::string traj_file_name) const {
  std::ofstream vmdf(file_name);
  vmdf << "display projection Orthographic" << std::endl
    << "color Display Background white" << std::endl
    << "axes location Off" << std::endl;
  vmdf << "topo readvarxyz " << trim("/", traj_file_name) << std::endl;
  vmdf << "mol modstyle 0 0 VDW 1.0000000 120.000000" << std::endl;
  double radius, distance;
  int center_index;
  for (int site_type = 0; site_type < config.num_site_types(); ++site_type) {
    vmdf << "set sel [atomselect top \"name " << site_type << "\"]" << std::endl;
    get_params(config, site_type, &radius, &distance, &center_index);
    vmdf << "$sel set radius " << radius << std::endl;
  }
}

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
  FEASST_CHECK_ALL_USED(args);
}

void FileXYZ::load(const std::string file_name, Configuration * config) const {
  std::ifstream xyz_file(file_name);
  ASSERT(xyz_file, "file_name: " << file_name << " is empty");
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
  TRACE("coord " << feasst::feasst_str(coord) << " pos " << position.str());
  config->set_side_lengths(position);
  // ASSERT(config->num_sites() == num_sites, "site mismatch");
  if (config->num_sites() == 0) {
    ASSERT(config->num_particle_types() > 0, "try adding particle type " <<
      "before loading the configuration");
    const int site_per_part = config->particle_types().particle(0).num_sites();
    ASSERT(num_sites % site_per_part == 0, "assumed particle type to load " <<
      "xyz file is incompatible with number of sites");
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
  if (euler_) {
    config->update_positions(coords, eulers);
  } else {
    config->update_positions(coords);
  }
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
  (*file.get()) << config.group_selects()[gindex].num_sites() << std::endl
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
