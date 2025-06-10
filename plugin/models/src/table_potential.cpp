#include <cmath>
#include <string>
#include <fstream>
#include "utils/include/arguments.h"
#include "utils/include/utils.h"
#include "utils/include/serialize.h"
#include "utils/include/io.h"
#include "math/include/constants.h"
#include "math/include/utils_math.h"
#include "configuration/include/model_params.h"
#include "configuration/include/domain.h"
#include "configuration/include/configuration.h"
#include "system/include/model_two_body.h"
#include "models/include/table_potential.h"

namespace feasst {

TablePotential::TablePotential(argtype * args) {
  class_name_ = "TablePotential";
  const std::string table_file = str("table_file", args, "-1");
  if (table_file != "-1") {
    read_table_(table_file);
  }
}
TablePotential::TablePotential(argtype args) : TablePotential(&args) {
  feasst_check_all_used(args);
}

void TablePotential::read_inner_(const int itype, const int jtype, const std::string& descript, const double value, std::ifstream& file) {
  ASSERT(descript == "inner", "format error: " << descript);
  const double inner = value;
  DEBUG("inner " << inner);
  ASSERT(inner >= 0, "inner: " << inner << " must be >= 0");
  inner_[itype][jtype] = inner;
  inner_g_[itype][jtype] = std::pow(inner, gamma_[itype][jtype]);
  if (inner_g_[itype][jtype] > 100) {
    WARN("inner(" << inner << ")^gamma(" << gamma_[itype][jtype] << ") is "
      << " large. Consider setting inner or gamma to a different value.");
  }
}

void TablePotential::read_table_(const std::string file_name) {
  DEBUG("file_name " << file_name);
  std::ifstream file(file_name);
  ASSERT(file.good(), "cannot find " << file_name);
  if (file.eof()) {
    return;
  }
  std::string line, descript;
  double double_val;
  int int_val;
  std::string str_val;
  file >> descript >> int_val;
  ASSERT(descript == "site_types", "format error: " << descript);
  const int num_sites = int_val;
  DEBUG("num_sites " << num_sites);

  // size arrays
  site_type_names_.resize(num_sites);
  resize(num_sites, num_sites, &inner_);
  resize(num_sites, num_sites, &inner_g_);
  resize(num_sites, num_sites, &cutoff_g_);
  resize(num_sites, num_sites, &energy_table_);
  resize(num_sites, num_sites, &gamma_);
  for (int type = 0; type < num_sites; ++type) {
    file >> str_val;
    DEBUG("site name " << str_val);
    site_type_names_[type] = str_val;
  }

  for (int itype = 0; itype < num_sites; ++itype) {
    for (int jtype = itype; jtype < num_sites; ++jtype) {
      // read optional "gamma" (default -2) and required "inner"
      file >> descript >> double_val;
      if (descript == "gamma") {
        const double gamma = double_val;
        DEBUG("gamma " << gamma);
        gamma_[itype][jtype] = gamma;

        // read inner
        file >> descript >> double_val;
        read_inner_(itype, jtype, descript, double_val, file);
      } else if (descript == "inner") {
        // If no gamma provided, set the default value
        gamma_[itype][jtype] = -2;
        read_inner_(itype, jtype, descript, double_val, file);
      } else {
        FATAL("table format error. Expected either \"gamma\" or \"inner\""
          << " but instead read: " << descript);
      }

      // read num_z
      file >> descript >> int_val;
      ASSERT(descript == "num_z", "format error: " << descript);
      const int num_z = int_val;
      DEBUG("num_z " << num_z);
      ASSERT(num_z > 0, "num_z: " << num_z << " must be > 0");

      // read energies
      Table1D * en = &energy_table_[itype][jtype];
      *en = Table1D({{"num", str(num_z)}});
      for (int z = 0; z < num_z; ++z) {
        file >> double_val;
        en->set_data(z, double_val);
      }
    }
  }

  // check for eof
  ASSERT(!file.eof(), "improper table file: " << file_name);
  file >> double_val;
  ASSERT(file.eof(), "improper table file: " << file_name);
}

void TablePotential::precompute(const Configuration& config) {
  Model::precompute(config);
  const ModelParams& existing = config.model_params();
  const ModelParam& cutoff = existing.select("cutoff");
  t2index_.assign(existing.size(), -1);
  site_types_.clear();
  for (const std::string& sname : site_type_names_) {
    site_types_.push_back(config.site_type_name_to_index(sname));
  }
  for (int t1 = 0; t1 < static_cast<int>(site_types_.size()); ++t1) {
    const int type1 = site_types_[t1];
    t2index_[type1] = t1;
  }
  for (int t1 = 0; t1 < static_cast<int>(site_types_.size()); ++t1) {
    const int type1 = site_types_[t1];
    for (int t2 = 0; t2 < static_cast<int>(site_types_.size()); ++t2) {
      const int type2 = site_types_[t2];
      const double rc = cutoff.mixed_value(type1, type2);
      cutoff_g_[t1][t2] = std::pow(rc, gamma_[t1][t2]);
    }
  }
}

double TablePotential:: energy(
    const double squared_distance,
    const int type1,
    const int type2,
    const ModelParams& model_params) {
  DEBUG("*** TablePotential ***");
  DEBUG("type1 " << type1 << " type2 " << type2);
  int ttype1 = type1;
  int ttype2 = type2;
  DEBUG("ttype1 " << ttype1 << " ttype2 " << ttype2);

  // enforce type1 <= type2 to avoid redundant tables.
  bool flip = false;
  if (ttype1 > ttype2) {
    flip = true;
  }
  if (flip) {
    feasst_swap(&ttype1, &ttype2);
  }
  DEBUG("flip " << flip);

  // convert site type to table type
  DEBUG("t2size " << t2index_.size());
  int tabtype1 = t2index_[ttype1];
  int tabtype2 = t2index_[ttype2];
  DEBUG("tabtype1 " << tabtype1 << " tabtype2 " << tabtype2);

  // Do not compute energy if either site is not represented.
  if (tabtype1 == -1 || tabtype2 == -1) {
    return 0.;
  }

  // check the inner cutoff.
  const double inner = inner_[tabtype1][tabtype2];
  DEBUG("squared_distance " << squared_distance);
  if (squared_distance < inner*inner) {
    TRACE("hard overlap");
    return NEAR_INFINITY;
  } else {
    const double gamma = gamma_[tabtype1][tabtype2];
    DEBUG("gamma " << gamma);
    const double rhg = inner_g_[tabtype1][tabtype2];
    DEBUG("rhg " << rhg);
    const double rcg = cutoff_g_[tabtype1][tabtype2];
    DEBUG("rcg " << rcg);
    double rg;
    if (gamma == -2) {
      rg = 1./squared_distance;
    } else {
      rg = std::pow(squared_distance, gamma/2.);
    }
    DEBUG("rg " << rg);
    const double z = (rg - rhg)/(rcg - rhg);
    DEBUG("z " << z);
    ASSERT(z >= 0 && z <= 1, "z: " << z);
    DEBUG("tab size " << energy_table_.size());
    DEBUG("tab size " << energy_table_[0].size());
    const double en = energy_table_[tabtype1][tabtype2].forward_difference_interpolation(z);
    DEBUG("en " << en);
    return en;
  }
}

FEASST_MAPPER(TablePotential,);

TablePotential::TablePotential(std::istream& istr) : ModelTwoBody(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version >= 8965 && version <= 8967, "unrecognized version: " << version);
  feasst_deserialize(&inner_, istr);
  feasst_deserialize(&inner_g_, istr);
  feasst_deserialize(&cutoff_g_, istr);
  feasst_deserialize(&site_types_, istr);
  if (version >= 8967) {
    feasst_deserialize(&site_type_names_, istr);
  }
  feasst_deserialize(&t2index_, istr);
  feasst_deserialize_fstobj(&energy_table_, istr);
  if (version >= 8966) {
    feasst_deserialize(&gamma_, istr);
  }
}

void TablePotential::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_model_(ostr);
  feasst_serialize_version(8967, ostr);
  feasst_serialize(inner_, ostr);
  feasst_serialize(inner_g_, ostr);
  feasst_serialize(cutoff_g_, ostr);
  feasst_serialize(site_types_, ostr);
  feasst_serialize(site_type_names_, ostr);
  feasst_serialize(t2index_, ostr);
  feasst_serialize_fstobj(energy_table_, ostr);
  feasst_serialize(gamma_, ostr);
}

}  // namespace feasst
