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

void TablePotential::read_inner_(const int itype, const int jtype, const std::string& descript, const double value) {
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
  std::getline(file, line);
  std::stringstream errmsg;
  errmsg << "unrecognized format of TablePotential::table_file (" << file_name
    << ") on line:" << line;
  std::vector<std::string> args = split(line, '=');
  int num_sites = -1, int_val;
  double double_val;
  if (args.size() == 1) {
    WARN("Deprecated table_file format with space-separated values.");
    args = split(line, ' ');
    std::string str_val;
    descript = args[0];
    //file >> descript >> int_val;
    ASSERT(descript == "site_types", "format error: " << descript);
    num_sites = str_to_int(args[1]);
    DEBUG("num_sites " << num_sites);
    ASSERT(static_cast<int>(args.size()) == num_sites + 2, errmsg.str());

    site_type_names_.resize(num_sites);
    for (int type = 0; type < num_sites; ++type) {
      //file >> str_val;
      str_val = args[2 + type];
      DEBUG("site name " << str_val);
      site_type_names_[type] = str_val;
    }
  } else if (args.size() == 2) {
    ASSERT(args[0] == "site_types", errmsg.str());
    site_type_names_ = split(args[1], ',');
    num_sites = static_cast<int>(site_type_names_.size());
  } else {
    FATAL(errmsg.str());
  }

  // size arrays
  resize(num_sites, num_sites, &inner_);
  resize(num_sites, num_sites, &inner_g_);
  resize(num_sites, num_sites, &cutoff_g_);
  resize(num_sites, num_sites, &energy_table_);
  resize(num_sites, num_sites, &gamma_);

  bool deprecated = false;
  for (int itype = 0; itype < num_sites; ++itype) {
    for (int jtype = itype; jtype < num_sites; ++jtype) {
      // If not provided, set the default value
      gamma_[itype][jtype] = -2;
      bool inner_set = false;

      DEBUG("itype:" << itype << " jtype:" << jtype);
      DEBUG("line: " << line);
      std::getline(file, line);
      args = split(line, '=');
      int num_z = -1;
      if (static_cast<int>(args.size()) == 1) {
        deprecated = true;
        WARN("Deprecated table_file format with space-separated values.");
        // read optional "gamma" (default -2) and required "inner"
        //file >> descript >> double_val;
        args = split(line, ' ');
        DEBUG("line:" << line);
        descript = args[0];
        double_val = str_to_double(args[1]);
        if (descript == "gamma") {
          const double gamma = double_val;
          DEBUG("gamma " << gamma);
          gamma_[itype][jtype] = gamma;

          // read inner
          file >> descript >> double_val;
          read_inner_(itype, jtype, descript, double_val);
        } else if (descript == "inner") {
          read_inner_(itype, jtype, descript, double_val);
        } else {
          FATAL("table format error. Expected either \"gamma\" or \"inner\""
            << " but instead read: " << descript);
        }

        // read num_z
        file >> descript >> int_val;
        ASSERT(descript == "num_z", "format error: " << descript);
        num_z = int_val;
        DEBUG("num_z " << num_z);
        ASSERT(num_z > 0, "num_z: " << num_z << " must be > 0");
      } else if (static_cast<int>(args.size()) == 2) {
        num_z = -1;
        bool read_data = false;
        while (!read_data) {
          DEBUG("args[0]:" << args[0]);
          if (args[0] == "gamma") {
            gamma_[itype][jtype] = str_to_double(args[1]);
            DEBUG("gamma " << gamma_[itype][jtype]);
          } else if (args[0] == "inner") {
            read_inner_(itype, jtype, "inner", str_to_double(args[1]));
            DEBUG("inner " << inner_[itype][jtype]);
            inner_set = true;
          } else if (args[0] == "num_z") {
            // No need to read num_z, but if set then check it
            num_z = str_to_int(args[1]);
            DEBUG("num_z " << num_z);
          } else {
            read_data = true;
          }
          if (!read_data) {
            std::getline(file, line);
            args = split(line, '=');
            DEBUG("line " << line);
          }
        }
        ASSERT(inner_set, "TablePotential requires inner argument not provided"
         << " by itype:" << itype << " jtype:" << jtype);
      } else {
        FATAL(errmsg.str());
      }
      // read energies
      std::vector<std::string> dat;
      if (!deprecated) {
        dat = split(line, ' ');
        if (num_z != -1) {
          ASSERT(num_z == static_cast<int>(dat.size()), "mismatch");
        } else {
          num_z = static_cast<int>(dat.size());
        }
      }
      Table1D * en = &energy_table_[itype][jtype];
      *en = Table1D({{"num", str(num_z)}});
      for (int z = 0; z < num_z; ++z) {
        if (deprecated) {
          file >> double_val;
        } else {
          double_val = str_to_double(dat[z]);
        }
        en->set_data(z, double_val);
      }
    }
  }

  // check for eof
  if (!file.eof()) {
    file >> double_val;
    ASSERT(file.eof(), "improper table file: " << file_name);
  }
}

void TablePotential::precompute(Configuration * config, ModelParams * params) {
  Model::precompute(config, params);
  const ModelParam& cutoff = params->select("cutoff");
  t2index_.assign(params->size(), -1);
  site_types_.clear();
  for (const std::string& sname : site_type_names_) {
    site_types_.push_back(config->site_type_name_to_index(sname));
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
      DEBUG("rc " << rc);
      DEBUG("gamma_[" << t1 << "][" << t2 << "]=" << gamma_[t1][t2]);
      cutoff_g_[t1][t2] = std::pow(rc, gamma_[t1][t2]);
      DEBUG("cutoff_g_[" << t1 << "][" << t2 << "]=" << cutoff_g_[t1][t2]);
    }
  }
}

double TablePotential:: energy(
    const double squared_distance,
    const int type1,
    const int type2,
    const ModelParams& model_params) {
  TRACE("*** TablePotential ***");
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
  const double inner = inner_[tabtype1][tabtype2];
  TRACE("squared_distance " << squared_distance);
  if (squared_distance < inner*inner) {
    TRACE("hard overlap");
    return NEAR_INFINITY;
  } else {
    const double gamma = gamma_[tabtype1][tabtype2];
    TRACE("gamma " << gamma);
    const double rhg = inner_g_[tabtype1][tabtype2];
    TRACE("rhg " << rhg);
    const double rcg = cutoff_g_[tabtype1][tabtype2];
    TRACE("rcg " << rcg);
    double rg;
    if (gamma == -2) {
      rg = 1./squared_distance;
    } else {
      rg = std::pow(squared_distance, gamma/2.);
    }
    TRACE("rg " << rg);
    const double z = (rg - rhg)/(rcg - rhg);
    TRACE("z " << z);
    ASSERT(z >= 0 && z <= 1, "z: " << z);
    TRACE("tab size " << energy_table_.size());
    TRACE("tab size " << energy_table_[0].size());
    const double en = energy_table_[tabtype1][tabtype2].forward_difference_interpolation(z);
    TRACE("en " << en);
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
