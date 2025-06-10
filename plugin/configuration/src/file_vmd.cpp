#include "utils/include/serialize.h"
#include "utils/include/io.h"
#include "utils/include/arguments.h"
#include "utils/include/debug.h"
#include "configuration/include/model_params.h"
#include "configuration/include/domain.h"
#include "configuration/include/file_vmd.h"

namespace feasst {

FileVMD::FileVMD(argtype * args) {
  min_sigma_ = dble("min_sigma", args, 0.1);
}
FileVMD::FileVMD(argtype args) : FileVMD(&args) {
  feasst_check_all_used(args);
}

void FileVMD::serialize(std::ostream& ostr) const {
  feasst_serialize_version(1, ostr);
  feasst_serialize(min_sigma_, ostr);
}

FileVMD::FileVMD(std::istream& istr) {
  feasst_deserialize_version(istr);
  feasst_deserialize(&min_sigma_, istr);
}

void FileVMD::get_params(const Configuration& config,
    const int site_type,
    double * radius,
    double * distance,
    int * center_index) const {
  double sigma = config.model_params().select("sigma").value(site_type);
  if (sigma < min_sigma_) {
    sigma = min_sigma_;
  }
  *radius = 0.5*sigma;
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
    vmdf << "set sel [atomselect top \"name " << site_type << "\"]"
         << std::endl;
    get_params(config, site_type, &radius, &distance, &center_index);
    vmdf << "$sel set radius " << radius << std::endl;
  }
  vmdf << "pbc set {";
  for (const double length : config.domain().side_lengths().coord()) {
    vmdf << length << " ";
  }
  vmdf << "} -all\npbc box -center origin -color blue";
}

}  // namespace feasst
