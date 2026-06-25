#include <iostream>
#include <cmath>
#include <string>
#include <fstream>
#include "utils/include/arguments.h"
#include "utils/include/utils.h"
#include "utils/include/serialize.h"
#include "utils/include/io.h"
#include "math/include/constants.h"
#include "math/include/utils_math.h"
#include "math/include/recursive_table.h"
#include "math/include/table.h"
#include "configuration/include/model_params.h"
#include "configuration/include/domain.h"
#include "configuration/include/configuration.h"
#include "system/include/model_two_body.h"
#include "aniso/include/model_recursive_table.h"

namespace feasst {

ModelRecursiveTable::ModelRecursiveTable(argtype * args) {
  class_name_ = "ModelRecursiveTable";
  input_file_ = str("input_file", args, "");
}
ModelRecursiveTable::ModelRecursiveTable(argtype args) : ModelRecursiveTable(&args) {
  feasst_check_all_used(args);
}
ModelRecursiveTable::~ModelRecursiveTable() {}

void files_to_types_recursive_(const std::vector<std::string>& files, const Configuration& config, std::vector<int> * t2index, std::vector<std::string> * site_type_names, std::vector<int> * site_types) {
  const ModelParams& existing = config.model_params();
  t2index->assign(existing.size(), -1);
  for (const std::string& filename : files) {
    const std::string splt = split(filename, ':')[0];
    const std::vector<std::string> types = split(splt, '_');
    ASSERT(static_cast<int>(types.size()) == 2,
      "input(" << splt << ") requires \"type1_type2:file\" format");
    for (const std::string& typ : types) {
      if (!find_in_list(typ, *site_type_names)) {
        site_type_names->push_back(typ);
      }
    }
  }
  for (const std::string& sname : *site_type_names) {
    site_types->push_back(config.site_type_name_to_index(sname));
  }
  for (int t1 = 0; t1 < static_cast<int>(site_types->size()); ++t1) {
    const int type1 = (*site_types)[t1];
    (*t2index)[type1] = t1;
  }
}

std::string filename_to_idx_recursive_(const std::string& filename, const Configuration& config, const std::vector<int>& t2index, int * idx1, int * idx2) {
  const std::vector<std::string> input = split(filename, ':');
  const std::vector<std::string> types = split(input[0], '_');
  ASSERT(static_cast<int>(types.size()) == 2,
    "input(" << input[0] << ") requires \"type1_type2:file\" format");
  const std::string inpfile = input[1];
  *idx1 = t2index[config.site_type_name_to_index(types[0])];
  *idx2 = t2index[config.site_type_name_to_index(types[1])];
  if (*idx1 > *idx2) feasst_swap(idx1, idx2);
  DEBUG("idx1:" << *idx1);
  DEBUG("idx2:" << *idx2);
  std::ifstream file(inpfile);
  std::string line;
  std::getline(file, line);
  ASSERT(!line.empty(), "file: " << inpfile << " is empty");
  return line;
}

void ModelRecursiveTable::precompute(Configuration * config, ModelParams * params) {
  if (static_cast<int>(energy_.size()) > 0) {
    return;
  }
  std::vector<std::string> site_type_names;
  ASSERT(!input_file_.empty(), "Error");
  std::vector<std::string> files = split(input_file_, ',');
  std::vector<int> site_types;
  files_to_types_recursive_(files, *config, &t2index_, &site_type_names, &site_types);

  const int num_sites = static_cast<int>(t2index_.size());
  resize(num_sites, num_sites, &lower_);
  resize(num_sites, num_sites, &upper_);
  resize(num_sites, num_sites, &energy_);
  DEBUG("t2index:" << feasst_str(t2index_));
  for (const std::string& filename : files) {
    int idx1, idx2;
    std::stringstream ss;
    ss << filename_to_idx_recursive_(filename, *config, t2index_, &idx1, &idx2);
    ModelRecursiveTable pot(ss);
    lower_[idx1][idx2] = pot.lower_[0][0];
    upper_[idx1][idx2] = pot.upper_[0][0];
    const double cutoff = pot.upper_[0][0];
    DEBUG("upper " << pot.upper_[0][0]);
    DEBUG("cutoff " << cutoff);
    const int type1 = site_types[idx1];
    const int type2 = site_types[idx2];
    params->set("cutoff", type1, type2, cutoff);
    params->set("cutoff", type2, type1, cutoff);
    std::cout << "# cutoff for " << type1 << "-" << type2 << " site types: " << cutoff << std::endl;
    std::cout << "# cutoff for " << type2 << "-" << type1 << " site types: " << cutoff << std::endl;
    //const double cut_conf = config->model_params().select("cutoff").mixed_value(type1, type2);
    //ASSERT(std::abs(cut_conf - cutoff) < 1e-7, "table cutoff: " << MAX_PRECISION << cutoff << " does not match configuration:" << cut_conf);
    energy_[idx1][idx2] = pot.energy_[0][0];
  }
}

double ModelRecursiveTable:: energy(
    const double squared_distance,
    const int type1,
    const int type2,
    const ModelParams& model_params) {
  TRACE("*** ModelRecursiveTable ***");
  TRACE("type1 " << type1 << " type2 " << type2);
  int ttype1 = type1;
  int ttype2 = type2;
  TRACE("ttype1 " << ttype1 << " ttype2 " << ttype2);

  // enforce type1 <= type2 to avoid redundant tables.
  bool flip = false;
  if (ttype1 > ttype2) {
    flip = true;
  }
  if (flip) {
    feasst_swap(&ttype1, &ttype2);
  }
  TRACE("flip " << flip);

  // convert site type to table type
  TRACE("t2size " << t2index_.size());
  int tabtype1 = t2index_[ttype1];
  int tabtype2 = t2index_[ttype2];
  TRACE("tabtype1 " << tabtype1 << " tabtype2 " << tabtype2);

  // Do not compute energy if either site is not represented.
  if (tabtype1 == -1 || tabtype2 == -1) {
    return 0.;
  }

  // check the inner cutoff.
  const double lwr = lower_[tabtype1][tabtype2];
  TRACE("lwr " << lwr);
  TRACE("squared_distance " << squared_distance);
  const double dist = std::sqrt(squared_distance);
  if (dist < lwr) {
    TRACE("hard overlap");
    return NEAR_INFINITY;
  } else {
    const double upr = upper_[tabtype1][tabtype2];
    TRACE("upr " << upr);
    const double z = (dist - lwr)/(upr - lwr);
    TRACE("z " << z);
    ASSERT(z >= 0 && z <= 1, "z: " << z);
    TRACE("tab size " << energy_.size());
    TRACE("tab size " << energy_[0].size());
    const double en = energy_[tabtype1][tabtype2].linear_interpolation(z);
    TRACE("en " << en);
    return en;
  }
}

FEASST_MAPPER(ModelRecursiveTable,);

ModelRecursiveTable::ModelRecursiveTable(std::istream& istr) : ModelTwoBody(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version >= 6765 && version <= 6765, "unrecognized version: " << version);
  feasst_deserialize(&lower_, istr);
  feasst_deserialize(&upper_, istr);
  feasst_deserialize(&t2index_, istr);
  feasst_deserialize_fstobj(&energy_, istr);
}

void ModelRecursiveTable::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_model_(ostr);
  feasst_serialize_version(6765, ostr);
  feasst_serialize(lower_, ostr);
  feasst_serialize(upper_, ostr);
  feasst_serialize(t2index_, ostr);
  feasst_serialize_fstobj(energy_, ostr);
}

}  // namespace feasst
