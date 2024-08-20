#include <fstream>
#include <sstream>
#include "utils/include/arguments.h"
#include "utils/include/utils.h"
#include "utils/include/serialize.h"
#include "utils/include/debug.h"
#include "configuration/include/domain.h"
#include "netcdf/include/file_netcdf.h"

namespace feasst {

FileNETCDF::FileNETCDF(argtype * args) : file_(file_parse_(args)) {
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
  euler_ = boolean("euler", args, false);
  ASSERT(!euler_, "euler isn't implemented yet");
  float_precision_ = boolean("float_precision", args, false);
}
FileNETCDF::FileNETCDF(argtype args) : FileNETCDF(&args) {
  feasst_check_all_used(args);
}

PrinterNETCDF::PrinterNETCDF(const bool euler,
    const bool float_precision,
    std::vector<int> * site_types,
    double * xyzs,
    float * xyzs_f) {
  euler_ = euler;
  float_precision_ = float_precision;
  site_types_ = site_types;
  xyzs_ = xyzs;
  xyzs_f_ = xyzs_f;
}

void PrinterNETCDF::work(const Site& site,
    const Configuration& config,
    const LoopDescriptor& data) {
  (*site_types_)[site_index_] = site.type();
  for (int dim = 0; dim < config.dimension(); ++dim) {
    if (xyzs_ == NULL) {
      *((xyzs_f_ + site_index_*config.dimension())+dim) = site.position().coord(dim);
    } else {
      *((xyzs_ + site_index_*config.dimension())+dim) = site.position().coord(dim);
    }
  }
  for (int dim = config.dimension(); dim < 3; ++dim) {
    if (xyzs_ == NULL) {
      *((xyzs_f_ + site_index_*config.dimension())+dim) = 0.;
    } else {
      *((xyzs_ + site_index_*config.dimension())+dim) = 0.;
    }
    //(*xyzs_)[site_index_][dim] = 0.;
  }
  if (euler_) {
//    const Euler& euler = site.euler();
//    (*file_.get()) << euler.phi() << " "
//                   << euler.theta() << " "
//                   << euler.psi();
  }
  ++site_index_;
}

void FileNETCDF::initialize(const Configuration& config) {
  const int num_dim = config.dimension();
  record_dim_ = file_->addDim("records");
  netCDF::NcDim nDim = file_->addDim("max_sites");
  netCDF::NcDim xDim = file_->addDim("dimensions", num_dim);
  xyz_dims_ = {record_dim_, nDim, xDim};
  domain_dims_ = {record_dim_, xDim};
  type_dims_ = {record_dim_, nDim};
  nc_domain_ = file_->addVar("domain_side_lengths", netCDF::ncDouble, domain_dims_);
  nc_tilt_xy_ = file_->addVar("domain_tilt_xy", netCDF::ncDouble, record_dim_);
  nc_tilt_xz_ = file_->addVar("domain_tilt_xz", netCDF::ncDouble, record_dim_);
  nc_tilt_yz_ = file_->addVar("domain_tilt_yz", netCDF::ncDouble, record_dim_);
  nc_site_types_ = file_->addVar("site_types", netCDF::ncInt, type_dims_);
  nc_num_sites_ = file_->addVar("num_sites", netCDF::ncInt, record_dim_);
  if (float_precision_) {
    nc_xyzs_ = file_->addVar("coordinates", netCDF::ncFloat, xyz_dims_);
  } else {
    nc_xyzs_ = file_->addVar("coordinates", netCDF::ncDouble, xyz_dims_);
  }
  num_config_ = 0;
}

void FileNETCDF::write(const Configuration& config) {
  const int gindex = gindex_(config);
  const int num_dim = config.dimension();
  const int num_site = config.num_sites(gindex);

  // add box length
  auto sls = config.domain().side_lengths().coord();
  double * side_lengths = &sls[0];
  std::vector<size_t> start = {static_cast<size_t>(num_config_), 0};
  std::vector<size_t> count = {1, static_cast<size_t>(num_dim)};
  //INFO("start " << feasst_str(start));
  nc_domain_.putVar(start, count, side_lengths);

  // add box tilts
  start = {static_cast<size_t>(num_config_)};
  count = {1};
  double xy = config.domain().xy();
  nc_tilt_xy_.putVar(start, count, &xy);
  double xz = config.domain().xz();
  nc_tilt_xz_.putVar(start, count, &xz);
  double yz = config.domain().yz();
  nc_tilt_yz_.putVar(start, count, &yz);

  // add number of sites
  nc_num_sites_.putVar(start, count, &num_site);

  // add site types
  std::vector<int> site_types(num_site, 0);

  // add coordinates
  double xyzs[num_site][num_dim];
  float xyzs_f[num_site][num_dim];
  if (float_precision_) {
    PrinterNETCDF printer(euler_, float_precision_, &site_types, NULL, (float *)xyzs_f);
    VisitConfiguration().loop(config, &printer, gindex);
  } else {
    PrinterNETCDF printer(euler_, float_precision_, &site_types, (double *)xyzs, NULL);
    VisitConfiguration().loop(config, &printer, gindex);
  }
  //INFO(feasst_str(site_types));
  int * site_types_p = &site_types[0];
  start = {static_cast<size_t>(num_config_), 0};
  count = {1, static_cast<size_t>(num_site)};
  nc_site_types_.putVar(start, count, site_types_p);
  start = {static_cast<size_t>(num_config_), 0, 0};
  count = {1, static_cast<size_t>(num_site), static_cast<size_t>(num_dim)};
  if (float_precision_) {
    nc_xyzs_.putVar(start, count, &xyzs_f);
  } else {
    nc_xyzs_.putVar(start, count, &xyzs);
  }
  ++num_config_;
}

void FileNETCDF::load_frame(std::ifstream& xyz_file, Configuration * config) const {
}

void FileNETCDF::load(const std::string file_name, Configuration * config) const {
//  std::ifstream xyz_file(file_name);
//  ASSERT(xyz_file, "file_name: " << file_name << " is empty");
//  load_frame(xyz_file, config);
}

void FileNETCDF::serialize(std::ostream& ostr) const {
  feasst_serialize_version(8856, ostr);
  feasst_serialize(group_index_, ostr);
  feasst_serialize(group_, ostr);
  feasst_serialize(euler_, ostr);
  feasst_serialize(float_precision_, ostr);
  feasst_serialize(num_config_, ostr);
}

FileNETCDF::FileNETCDF(std::istream& istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 8856, "version mismatch: " << version);
  feasst_deserialize(&group_index_, istr);
  feasst_deserialize(&group_, istr);
  feasst_deserialize(&euler_, istr);
  feasst_deserialize(&float_precision_, istr);
  feasst_deserialize(&num_config_, istr);
}

}  // namespace feasst
