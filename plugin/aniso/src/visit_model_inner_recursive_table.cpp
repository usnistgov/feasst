#include <cmath>  // isnan, pow
#include <string>
#include <fstream>
#include "utils/include/utils.h"  // resize
#include "utils/include/serialize.h"
#include "utils/include/arguments.h"
#include "math/include/utils_math.h"
#include "math/include/constants.h"
#include "configuration/include/particle_factory.h"
#include "configuration/include/configuration.h"
#include "aniso/include/model_recursive_table.h"
#include "aniso/include/visit_model_inner_recursive_table.h"

namespace feasst {

VisitModelInnerRecursiveTable::VisitModelInnerRecursiveTable(argtype * args) : VisitModelInnerTable(args) {
  class_name_ = "VisitModelInnerRecursiveTable";
  input_file_ = str("input_file", args, "");
}
VisitModelInnerRecursiveTable::VisitModelInnerRecursiveTable(argtype args) : VisitModelInnerRecursiveTable(&args) {
  feasst_check_all_used(args);
}

void VisitModelInnerRecursiveTable::read_table(const std::string file_name,
    const bool ignore_energy,
    Configuration * config) {
  if (static_cast<int>(contact_.size()) > 0 ||
      static_cast<int>(contact2d_.size()) > 0) {
    return;
  }
  std::vector<std::string> site_type_names;
  ASSERT(!input_file_.empty(), "Error");
  std::vector<std::string> files = split(input_file_, ',');
  std::vector<int> site_types;
  files_to_types_recursive_(files, *config, &t2index_, &site_type_names, &site_types);

  const int num_sites = static_cast<int>(t2index_.size());
  //resize(num_sites, num_sites, &lower_);
  //resize(num_sites, num_sites, &upper_);
  resize(num_sites, num_sites, &contact_);
  resize(num_sites, num_sites, &contact2d_);
  DEBUG("t2index:" << feasst_str(t2index_));
  for (const std::string& filename : files) {
    int idx1, idx2;
    std::stringstream ss;
    ss << filename_to_idx_recursive_(filename, *config, t2index_, &idx1, &idx2);
    VisitModelInnerRecursiveTable vis(ss);
    //lower_[idx1][idx2] = vis.lower_[0][0];
    //upper_[idx1][idx2] = vis.upper_[0][0];
    //const double cutoff = vis.upper_[0][0];
    //DEBUG("upper " << vis.upper_[0][0]);
    //DEBUG("cutoff " << cutoff);
    //const int type1 = site_types[idx1];
    //const int type2 = site_types[idx2];
    //config->set_model_param("cutoff", type1, type2, cutoff);
    //config->set_model_param("cutoff", type2, type1, cutoff);
    //std::cout << "# cutoff for " << type1 << "-" << type2 << " site types: " << cutoff << std::endl;
    //std::cout << "# cutoff for " << type2 << "-" << type1 << " site types: " << cutoff << std::endl;
    //const double cut_conf = config->model_params().select("cutoff").mixed_value(type1, type2);
    //ASSERT(std::abs(cut_conf - cutoff) < 1e-7, "table cutoff: " << MAX_PRECISION << cutoff << " does not match configuration:" << cut_conf);
    if (vis.contact_.size() > 0) {
      contact_[idx1][idx2] = vis.contact_[0][0];
    } else {
      contact2d_[idx1][idx2] = vis.contact2d_[0][0];
    }
  }
}

double VisitModelInnerRecursiveTable::compute_aniso(const int type1, const int type2,
    const double squared_distance, const double s1, const double s2,
    const Configuration& config) const {
  TRACE("size1 " << contact2d_.size());
  TRACE("size2 " << contact2d_[0].size());
  const float inner = contact2d_[type1][type2].linear_interpolation(s1, s2);
  TRACE("inner " << inner << " " << inner*inner);
  TRACE("squared_distance " << squared_distance);
  double en = 0.;
  if (squared_distance < inner*inner) {
    en = NEAR_INFINITY;
    TRACE("hard overlap");
  } else if (ignore_energy_) {
    en = 0.;
  } else {
    FATAL("implement");
  }
  return en;
}

double VisitModelInnerRecursiveTable::compute_aniso(const int type1, const int type2,
    const double squared_distance, const double s1, const double s2,
    const double e1, const double e2, const double e3, const Configuration& config) const {
  TRACE("size1 " << contact_.size());
  TRACE("size2 " << contact_[0].size());
  const float inner = contact_[type1][type2].linear_interpolation(s1, s2, e1, e2, e3);
  TRACE("inner " << inner << " " << inner*inner);
  TRACE("squared_distance " << squared_distance);
  double en = 0.;
  if (squared_distance < inner*inner) {
    en = NEAR_INFINITY;
    TRACE("hard overlap");
  } else if (ignore_energy_) {
    en = 0.;
  } else {
    FATAL("implement");
//    const double delta = delta_[type1][type2];
//    const double outer = inner + delta;
//    TRACE("delta " << delta);
//    TRACE("outer " << outer);
//    if (squared_distance < outer*outer) {
//      const double gamma = gamma_[type1][type2];
//      TRACE("gamma " << gamma);
//      const std::vector<std::vector<std::shared_ptr<Table6D> > >& energyt = config.table6d();
//      if ((std::abs(gamma) < NEAR_ZERO)) {
//        en = -1;
//      } else if (is_energy_table(energyt)) {
//        const double smooth = smoothing_distance_[type1][type2];
//        const double rhg = std::pow(inner, gamma);
//        const double rcg = std::pow(outer - smooth, gamma);
//        const double rg = std::pow(squared_distance, 0.5*gamma);
//        double z = (rg - rhg)/(rcg - rhg);
//        if (z < 0 && z > -1e-6) {
//          z = 0.;
//        }
//        TRACE("z " << z);
//        if (z > 1.) {
//          en = energyt[type1][type2]->linear_interpolation(s1, s2, e1, e2, e3, 1.);
//          const double dx = outer - std::sqrt(squared_distance);
//          TRACE("dx " << dx);
//          if (dx > smooth && dx < smooth + 1e-5) {
//            en = 0.;
//          } else {
//            ASSERT(dx >= 0 && dx <= smooth, "dx: " << MAX_PRECISION << dx);
//            en *= dx/smooth;
//          }
//        } else {
//          ASSERT(z >= 0 && z <= 1, "z: " << MAX_PRECISION << z);
//          en = energyt[type1][type2]->linear_interpolation(s1, s2, e1, e2, e3, z);
//        }
//      }
//    }
  }
  return en;
}

FEASST_MAPPER(VisitModelInnerRecursiveTable,);

VisitModelInnerRecursiveTable::VisitModelInnerRecursiveTable(std::istream& istr) : VisitModelInnerTable(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 6855, "unrecognized version: " << version);
  feasst_deserialize(&input_file_, istr);
  feasst_deserialize_fstobj(&contact_, istr);
  feasst_deserialize_fstobj(&contact2d_, istr);
}

void VisitModelInnerRecursiveTable::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_visit_model_inner_table_(ostr);
  feasst_serialize_version(6855, ostr);
  feasst_serialize(input_file_, ostr);
  feasst_serialize_fstobj(contact_, ostr);
  feasst_serialize_fstobj(contact2d_, ostr);
}

}  // namespace feasst
