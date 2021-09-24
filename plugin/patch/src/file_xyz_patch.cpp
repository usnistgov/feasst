
#include <fstream>
#include <sstream>
#include "utils/include/serialize.h"
#include "utils/include/debug.h"
#include "math/include/utils_math.h"
#include "configuration/include/domain.h"
#include "patch/include/file_xyz_patch.h"

namespace feasst {

void FileVMDPatch::serialize(std::ostream& ostr) const {
  feasst_serialize_version(3284, ostr);
}

FileVMDPatch::FileVMDPatch(std::istream& istr) : FileVMD(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 3284, "version mismatch: " << version);
}

void FileVMDPatch::get_params(const Configuration& config,
    const int site_type,
    double * radius,
    double * distance,
    int *center_index) const {
  if (config.model_params().select("director").value(site_type) < 0.5) {
    const double sigma = config.model_params().sigma().value(site_type);
    *radius = 0.5*sigma;
    *distance = 0.;
    *center_index = -1;
  } else {
    const double pa = config.model_params().select("patch_angle").value(site_type);
    const double cpa = std::cos(degrees_to_radians(pa));
    //const double cpa = config.model_params().select("cos_patch_angle").value(site_type);
    DEBUG("cpa " << cpa);

    // find first bond of site of type to get sigma
    int ptype = 0;
    int prev_sts = 0;
    while (prev_sts + config.unique_type(ptype).num_sites() < site_type) {
      prev_sts += config.unique_type(ptype).num_sites();
      ++ptype;
    }
    DEBUG("ptype " << ptype);
    DEBUG("prev_sts " << prev_sts);
    int site_index = 0;
    while (config.particle_type(ptype).site(site_index).type() != site_type) {
      ++site_index;
    }
    DEBUG("site_index " << site_index);
    *center_index = config.particle_type(ptype).bond_neighbors(site_index)[0];
    DEBUG("center_index " << *center_index);
    const int center_type = config.particle_type(ptype).site(*center_index).type();
    DEBUG("center_type " << center_type);
    const double sigma = config.model_params().sigma().value(center_type);

    const double r1 = sigma/2, eps = r1/20.;
    DEBUG("r1 " << r1 << " eps " << eps);
    DEBUG("cpa " << cpa);
    *radius = (2*r1*(1-cpa)*(r1+eps)+eps*eps)/(2*eps+2*r1*(1-cpa));
    DEBUG("radius " << *radius);
    *distance = r1 + eps - *radius;
    DEBUG("distance " << *distance);
  }
}

FileXYZPatch::FileXYZPatch(argtype * args) {
  group_index_ = integer("group_index", args, 0);
  append_ = boolean("append", args, false);
}
FileXYZPatch::FileXYZPatch(argtype args) : FileXYZPatch(&args) {
  check_all_used(args);
}

void FileXYZPatch::load(const std::string file_name, Configuration * config) const {
  FATAL("not implemented");
}

PrinterXYZPatch::PrinterXYZPatch(std::shared_ptr<std::ofstream> file,
    const int num_places) : file_(file) {
  num_places_ = num_places;
}

void PrinterXYZPatch::work(const Site& site,
    const Configuration& config,
    const LoopDescriptor& data) {
  (*file_.get()) << site.type() << " ";
  (*file_.get()) << std::setprecision(num_places_);
  double radius, distance;
  int center_index;
  vmd_.get_params(config, site.type(), &radius, &distance, &center_index);
  DEBUG("center index " << center_index);
  if (center_index < 0) {
    for (int dim = 0; dim < config.dimension(); ++dim) {
      (*file_.get()) << site.position().coord(dim) << " ";
    }
  } else {
    Position center = config.select_particle(data.particle_index).site(center_index).position();
    Position diff = site.position();
    diff.subtract(center);
    diff.normalize();
    diff.multiply(distance);
    diff.add(center);
    for (int dim = 0; dim < config.dimension(); ++dim) {
      (*file_.get()) << diff.coord(dim) << " ";
    }
  }
  for (int dim = config.dimension(); dim < 3; ++dim) {
    (*file_.get()) << "0 ";
  }
  (*file_.get()) << std::endl;
}

void FileXYZPatch::write(const std::string file_name,
                    const Configuration& config,
                    const int num_places) const {
  auto file = std::make_shared<std::ofstream>();
  if (append_ == 0) {
    file->open(file_name);
  } else {
    file->open(file_name, std::ofstream::app);
  }
  const Domain& domain = config.domain();
  (*file.get()) << config.group_selects()[group_index_].num_sites() << std::endl
    << "-1 ";
  (*file.get()) << std::setprecision(num_places);
  for (int dim = 0; dim < domain.dimension(); ++dim) {
    (*file.get()) << domain.side_length(dim) << " ";
  }
  (*file.get()) << domain.xy() << " "
    << domain.xz() << " "
    << domain.yz() << " "
    << std::endl;
  PrinterXYZPatch printer(file);
  VisitConfiguration().loop(config, &printer, group_index_);
}

void FileXYZPatch::write_for_vmd(const std::string file_name,
    const Configuration& config) const {
  std::stringstream ss;
  ss << file_name << ".vmd";
  write(file_name, config);
  FileVMDPatch().write(ss.str(), config, file_name);
}

void FileXYZPatch::serialize(std::ostream& ostr) const {
  feasst_serialize_version(6908, ostr);
  feasst_serialize(group_index_, ostr);
  feasst_serialize(append_, ostr);
}

FileXYZPatch::FileXYZPatch(std::istream& istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 6908, "version mismatch: " << version);
  feasst_deserialize(&group_index_, istr);
  feasst_deserialize(&append_, istr);
}
}  // namespace feasst
