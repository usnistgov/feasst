#include <iostream>
#include <cmath>  // isnan, pow
#include <string>
#include <fstream>
#include "utils/include/utils.h"  // resize
#include "utils/include/serialize.h"
#include "utils/include/arguments.h"
#include "math/include/utils_math.h"
#include "math/include/constants.h"
#include "configuration/include/domain.h"
#include "configuration/include/particle_factory.h"
#include "configuration/include/configuration.h"
#include "configuration/include/model_params.h"
#include "aniso/include/model_recursive_table.h"
#include "aniso/include/recursive_table.h"

namespace feasst {

RecursiveTable::RecursiveTable(argtype * args) : VisitModelInnerTable(args) {
  class_name_ = "RecursiveTable";
  input_file_ = str("input_file", args, "");
}
RecursiveTable::RecursiveTable(argtype args) : RecursiveTable(&args) {
  feasst_check_all_used(args);
}

void RecursiveTable::read_table(const std::string file_name,
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
  resize(num_sites, num_sites, &contact1d_);
  resize(num_sites, num_sites, &cutoff_);
  resize(num_sites, num_sites, &cutoff2d_);
  resize(num_sites, num_sites, &cutoff1d_);
  resize(num_sites, num_sites, &energy_);
  resize(num_sites, num_sites, &energy3d_);
  resize(num_sites, num_sites, &energy2d_);
  DEBUG("t2index:" << feasst_str(t2index_));
  for (const std::string& filename : files) {
    int idx1, idx2;
    std::stringstream ss;
    ss << filename_to_idx_recursive_(filename, *config, t2index_, &idx1, &idx2);
    RecursiveTable vis(ss);
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

    // read the tables and obtain maximum contact or cutoff values
    DEBUG("vis.contact_.size():" << vis.contact_.size());
    DEBUG("vis.contact2d_.size():" << vis.contact2d_.size());
    DEBUG("vis.contact1d_.size():" << vis.contact1d_.size());
    double max_contact = -1;
    if (vis.contact_.size() > 0) {
      contact_[idx1][idx2] = vis.contact_[0][0];
      max_contact = contact_[idx1][idx2].maximum();
    } else if (vis.contact2d_.size() > 0) {
      contact2d_[idx1][idx2] = vis.contact2d_[0][0];
      max_contact = contact2d_[idx1][idx2].maximum();
    } else if (vis.contact1d_.size() > 0) {
      contact1d_[idx1][idx2] = vis.contact1d_[0][0];
      max_contact = contact1d_[idx1][idx2].maximum();
    }

    double max_cutoff = -1;
    DEBUG("vis.cutoff_.size():" << vis.cutoff_.size());
    DEBUG("vis.cutoff2d_.size():" << vis.cutoff2d_.size());
    DEBUG("vis.cutoff1d_.size():" << vis.cutoff1d_.size());
    if (vis.cutoff_.size() > 0) {
      cutoff_[idx1][idx2] = vis.cutoff_[0][0];
      max_cutoff = cutoff_[idx1][idx2].maximum();
    } else if (vis.cutoff2d_.size() > 0) {
      cutoff2d_[idx1][idx2] = vis.cutoff2d_[0][0];
      max_cutoff = cutoff2d_[idx1][idx2].maximum();
    } else if (vis.cutoff1d_.size() > 0) {
      DEBUG("here");
      cutoff1d_[idx1][idx2] = vis.cutoff1d_[0][0];
      max_cutoff = cutoff1d_[idx1][idx2].maximum();
      DEBUG("max_cutoff " << max_cutoff);
    }

    if (vis.energy_.size() > 0) {
      energy_[idx1][idx2] = vis.energy_[0][0];
    } else if (vis.energy3d_.size() > 0) {
      energy3d_[idx1][idx2] = vis.energy3d_[0][0];
    } else if (vis.energy2d_.size() > 0) {
      energy2d_[idx1][idx2] = vis.energy2d_[0][0];
    }

    // automatically set the cutoff
    double cutoff = max_cutoff;
    DEBUG("max_cutoff " << max_cutoff);
    DEBUG("ignore_energy: " << ignore_energy);
    if (max_cutoff == -1) {
      cutoff = max_contact;
      std::cout << "# RecursiveTable::ignore_energy=true set automatically." << std::endl;
      ignore_energy_ = true;
      //ASSERT(ignore_energy_, "energy was ignored, so using contact as cutoff.");
    } else {
      ASSERT(!ignore_energy_, "energy is not ignored, so should be using cutoff.");
    }
    ASSERT(cutoff > 0, "err");
    const int type1 = site_types[idx1];
    const int type2 = site_types[idx2];
    ASSERT(cutoff <= config->domain().min_side_length()/2., "the automatically "
      << "computed cutoff:" << cutoff << " is > half the minimum side length:"
      << config->domain().min_side_length()/2.);
    config->set_model_param("cutoff", type1, type2, cutoff);
    config->set_model_param("cutoff", type2, type1, cutoff);
    std::cout << "# RecursiveTable cutoff for " << type1 << "-" << type2 << " site types: " << cutoff << std::endl;
    std::cout << "# RecursiveTable cutoff for " << type2 << "-" << type1 << " site types: " << cutoff << std::endl;
  }
}

double RecursiveTable::compute_aniso(const int tabtype1, const int tabtype2,
    const double squared_distance, const double s1,
    const Configuration& config, const ModelParams& model_params,
    const int stype1, const int stype2) const {
  const double sigma = model_params.select(sigma_index()).mixed_values()[stype1][stype2];
  ASSERT(std::abs(sigma - 1) < 1e-8, sigma);
  const float inner = sigma*contact1d_[tabtype1][tabtype2].linear_interpolation(s1);
  TRACE("inner " << inner << " " << inner*inner);
  TRACE("squared_distance " << squared_distance);
  double en = 0.;
  if (squared_distance < inner*inner) {
    en = NEAR_INFINITY;
    TRACE("hard overlap");
  } else if (!ignore_energy_) {
    const float cutoff = sigma*cutoff1d_[tabtype1][tabtype2].linear_interpolation(s1);
    TRACE("cutoff " << cutoff);
    if (squared_distance < cutoff*cutoff) {
      const double epsilon = model_params.select(epsilon_index()).mixed_values()[stype1][stype2];
      float z = (std::sqrt(squared_distance)-inner)/(cutoff - inner);
      if (z < 0. && z > -1e-6) {
        z = 0.;
      }
      TRACE("z " << z);
      if (z <= 1.) {
        en = epsilon*energy2d_[tabtype1][tabtype2].linear_interpolation(s1, z);
      }
    }
  }
  return en;
}

double RecursiveTable::compute_aniso(const int tabtype1, const int tabtype2,
    const double squared_distance, const double s1, const double s2,
    const Configuration& config, const ModelParams& model_params,
    const int stype1, const int stype2) const {
  TRACE("size1 " << contact2d_.size());
  TRACE("size2 " << contact2d_[0].size());
  TRACE("s1 " << s1 << " s2 " << s2);
  const double sigma = model_params.select(sigma_index()).mixed_values()[stype1][stype2];
  const float inner = sigma*contact2d_[tabtype1][tabtype2].linear_interpolation(s1, s2);
  TRACE("inner " << inner << " " << inner*inner);
  TRACE("squared_distance " << squared_distance);
  double en = 0.;
  if (squared_distance < inner*inner) {
    en = NEAR_INFINITY;
    TRACE("hard overlap");
  } else if (!ignore_energy_) {
    const float cutoff = sigma*cutoff2d_[tabtype1][tabtype2].linear_interpolation(s1, s2);
    if (squared_distance < cutoff*cutoff) {
      const double epsilon = model_params.select(epsilon_index()).mixed_values()[stype1][stype2];
      float z = (std::sqrt(squared_distance)-inner)/(cutoff - inner);
      if (z < 0. && z > -1e-6) {
        z = 0.;
      }
      TRACE("z " << z);
      if (z <= 1.) {
        en = epsilon*energy3d_[tabtype1][tabtype2].linear_interpolation(s1, s2, z);
      }
    }
  }
  return en;
}

double RecursiveTable::compute_aniso(const int tabtype1, const int tabtype2,
    const double squared_distance, const double s1, const double s2,
    const double e1, const double e2, const double e3,
    const Configuration& config, const ModelParams& model_params,
    const int stype1, const int stype2) const {
  TRACE("size1 " << contact_.size());
  TRACE("size2 " << contact_[0].size());
  const double sigma = model_params.select(sigma_index()).mixed_values()[stype1][stype2];
  const float inner = sigma*contact_[tabtype1][tabtype2].linear_interpolation(s1, s2, e1, e2, e3);
  TRACE("inner " << inner << " " << inner*inner);
  TRACE("squared_distance " << squared_distance);
  double en = 0.;
  if (squared_distance < inner*inner) {
    en = NEAR_INFINITY;
    TRACE("hard overlap");
  } else if (!ignore_energy_) {
    const float cutoff = sigma*cutoff_[tabtype1][tabtype2].linear_interpolation(s1, s2, e1, e2, e3);
    if (squared_distance < cutoff*cutoff) {
      const double epsilon = model_params.select(epsilon_index()).mixed_values()[stype1][stype2];
      float z = (std::sqrt(squared_distance)-inner)/(cutoff - inner);
      if (z < 0. && z > -1e-6) {
        z = 0.;
      }
      TRACE("z " << z);
      if (z <= 1.) {
        en = epsilon*energy_[tabtype1][tabtype2].linear_interpolation(s1, s2, e1, e2, e3, z);
      }
    }
  }
  return en;
}

FEASST_MAPPER(RecursiveTable,);

RecursiveTable::RecursiveTable(std::istream& istr) : VisitModelInnerTable(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 6855, "unrecognized version: " << version);
  feasst_deserialize(&input_file_, istr);
  feasst_deserialize_fstobj(&contact_, istr);
  feasst_deserialize_fstobj(&contact2d_, istr);
  feasst_deserialize_fstobj(&cutoff_, istr);
  feasst_deserialize_fstobj(&cutoff2d_, istr);
  feasst_deserialize_fstobj(&energy_, istr);
  feasst_deserialize_fstobj(&energy3d_, istr);
  feasst_deserialize_fstobj(&contact1d_, istr);
  feasst_deserialize_fstobj(&cutoff1d_, istr);
  feasst_deserialize_fstobj(&energy2d_, istr);
}

void RecursiveTable::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_visit_model_inner_table_(ostr);
  feasst_serialize_version(6855, ostr);
  feasst_serialize(input_file_, ostr);
  feasst_serialize_fstobj(contact_, ostr);
  feasst_serialize_fstobj(contact2d_, ostr);
  feasst_serialize_fstobj(cutoff_, ostr);
  feasst_serialize_fstobj(cutoff2d_, ostr);
  feasst_serialize_fstobj(energy_, ostr);
  feasst_serialize_fstobj(energy3d_, ostr);
  feasst_serialize_fstobj(contact1d_, ostr);
  feasst_serialize_fstobj(cutoff1d_, ostr);
  feasst_serialize_fstobj(energy2d_, ostr);
}

}  // namespace feasst
