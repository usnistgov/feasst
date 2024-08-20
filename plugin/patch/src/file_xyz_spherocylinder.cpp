#include <fstream>
#include "utils/include/arguments.h"
#include "utils/include/serialize.h"
#include "utils/include/debug.h"
#include "math/include/utils_math.h"
#include "configuration/include/select.h"
#include "configuration/include/particle_factory.h"
#include "configuration/include/model_params.h"
#include "configuration/include/domain.h"
#include "patch/include/file_xyz_spherocylinder.h"

namespace feasst {

void FileVMDSpherocylinder::serialize(std::ostream& ostr) const {
  FileVMD::serialize(ostr);
  feasst_serialize_version(5740, ostr);
}

FileVMDSpherocylinder::FileVMDSpherocylinder(std::istream& istr) : FileVMD(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 5740, "version mismatch: " << version);
}

void FileVMDSpherocylinder::get_params(const Configuration& config,
    const int site_type,
    double * radius,
    double * length,
    int *center_index) const {
  if (config.model_params().select("director").value(site_type) < 0.5) {
    const double sigma = config.model_params().select("sigma").value(site_type);
    *radius = 0.5*sigma;
    *length = 0.;
    *center_index = -1;
  } else {
    *length = config.model_params().select("spherocylinder_length").value(site_type);
    DEBUG("length " << *length);

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
    const double sigma = config.model_params().select("sigma").value(center_type);
    *radius = 0.5*sigma;
    DEBUG("radius " << *radius);
  }
}

FileXYZSpherocylinder::FileXYZSpherocylinder(argtype * args) {
  expand_factor_ = integer("expand_factor", args, 10);
  group_index_ = integer("group_index", args, 0);
  append_ = boolean("append", args, false);
}
FileXYZSpherocylinder::FileXYZSpherocylinder(argtype args) : FileXYZSpherocylinder(&args) {
  feasst_check_all_used(args);
}

void FileXYZSpherocylinder::load(const std::string file_name, Configuration * config) const {
  FATAL("not implemented");
}

PrinterXYZSpherocylinder::PrinterXYZSpherocylinder(std::shared_ptr<std::ofstream> file,
    const int num_places, const int expand_factor) : file_(file) {
  num_places_ = num_places;
  expand_factor_ = expand_factor;
}

void PrinterXYZSpherocylinder::work(const Site& site,
    const Configuration& config,
    const LoopDescriptor& data) {
  double radius, length;
  int center_index;
  vmd_.get_params(config, site.type(), &radius, &length, &center_index);
  DEBUG("center index " << center_index);
  if (center_index >= 0) {
    Position center = config.select_particle(data.particle_index).site(center_index).position();
    //Position end2 = end1;
    const double dz = 1./(expand_factor_-1);
    for (int iz = 0; iz < expand_factor_; ++iz) {
      const double z = -0.5+iz*dz;
      //for (double z = -0.5; z <= 0.5-1e-8; z += dz) {
      Position end1 = site.position();
      end1.subtract(center);
      end1.normalize();
      end1.multiply(z*length);
      //end2.multiply(-0.5*length);
      end1.add(center);
      //end2.add(center);
      (*file_.get()) << site.type() << " " << std::setprecision(num_places_);
      for (int dim = 0; dim < config.dimension(); ++dim) {
        (*file_.get()) << end1.coord(dim) << " ";
      }
      for (int dim = config.dimension(); dim < 3; ++dim) {
        (*file_.get()) << "0 ";
      }
      (*file_.get()) << std::endl;
      //<< site.type() << " ";
//      for (int dim = 0; dim < config.dimension(); ++dim) {
//        (*file_.get()) << end2.coord(dim) << " ";
//      }
//      for (int dim = config.dimension(); dim < 3; ++dim) {
//        (*file_.get()) << "0 ";
//      }
//      (*file_.get()) << std::endl;
    }
  }
}

void FileXYZSpherocylinder::write(const std::string file_name,
                    const Configuration& config,
                    const int num_places) const {
  auto file = std::make_shared<std::ofstream>();
  if (append_ == 0) {
    file->open(file_name);
  } else {
    file->open(file_name, std::ofstream::app);
  }
  const Domain& domain = config.domain();
  (*file.get()) << config.group_select(group_index_).num_sites()*expand_factor_/2 << std::endl
    << "-1 ";
  (*file.get()) << std::setprecision(num_places);
  for (int dim = 0; dim < domain.dimension(); ++dim) {
    (*file.get()) << domain.side_length(dim) << " ";
  }
  (*file.get()) << domain.xy() << " "
    << domain.xz() << " "
    << domain.yz() << " "
    << std::endl;
  PrinterXYZSpherocylinder printer(file, 8, expand_factor_);
  VisitConfiguration().loop(config, &printer, group_index_);
}

void FileXYZSpherocylinder::write_for_vmd(const std::string file_name,
    const Configuration& config) const {
  std::stringstream ss;
  ss << file_name << ".vmd";
  write(file_name, config);
  FileVMDSpherocylinder().write(ss.str(), config, file_name);
}

void FileXYZSpherocylinder::serialize(std::ostream& ostr) const {
  feasst_serialize_version(2034, ostr);
  feasst_serialize(expand_factor_, ostr);
  feasst_serialize(group_index_, ostr);
  feasst_serialize(append_, ostr);
}

FileXYZSpherocylinder::FileXYZSpherocylinder(std::istream& istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 2034, "version mismatch: " << version);
  feasst_deserialize(&expand_factor_, istr);
  feasst_deserialize(&group_index_, istr);
  feasst_deserialize(&append_, istr);
}
}  // namespace feasst
