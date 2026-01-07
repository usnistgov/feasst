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
#include "models/include/recursive_table_potential.h"

namespace feasst {

RecursiveTablePotential::RecursiveTablePotential(argtype * args) {
  class_name_ = "RecursiveTablePotential";
  mayer_training_file_ = str("mayer_training_file", args, "");
  hard_limit_u_ = dble("hard_limit_u", args, 100);
  size_ = dble("size", args, 16);
  beta_ = dble("beta", args, 1.);
  min_criteria_ = dble("min_criteria", args, 0.03);
}
RecursiveTablePotential::RecursiveTablePotential(argtype args) : RecursiveTablePotential(&args) {
  feasst_check_all_used(args);
}
RecursiveTablePotential::~RecursiveTablePotential() {}

void RecursiveTablePotential::build_table_(const double lower, const double upper, const std::vector<std::vector<double> >& data, RecursiveTable1D * tab) const {
  for (int i = 0; i < tab->num(); ++i) {
    // find the two data points closest to the value of interest
    const double z = static_cast<double>(i)/(tab->num() - 1);
    const double r2 = z2r_(z, lower, upper);
    DEBUG("r2:" << r2);
    bool found = false;
    int index = -1;
    while (!found) {
      ++index;
      if (data[index][0] > r2) {
        DEBUG("found " << index);
        found = true;
        ASSERT(index != 0, "er");
      }
      ASSERT(index < static_cast<int>(data.size()), "cannot find");
    }
    // linearly interpolate between the data points to set the data table
    const double x0 = data[index-1][0];
    const double x1 = data[index][0];
    const double y0 = data[index-1][1];
    const double y1 = data[index][1];
    DEBUG("x0:" << x0 << " x1:" << x1 << " y0:" << y0 << " y1:" << y1);
    const double y = y0 + (y1-y0)*(1.-(x1 - r2)/(x1 - x0));
    DEBUG("y: " << y);
    tab->set_data(i, y);
  }
}

int RecursiveTablePotential::analyze_table_(const double lower, const double upper, const std::vector<std::vector<double> >& data, const Table1D& table, const double beta, const std::string& filename, double * max_criteria) {
  // compare the table with the data
  // find the global max of |e^-bU - e^-bU_table| to insert a recurive table at that point
  //for (const auto& dat : data) {
  //int max = -1;
  double max_z = -1;
  *max_criteria = -1;
  std::ofstream file(filename);
  for (int idat = 0; idat < static_cast<int>(data.size()); ++idat) {
    const std::vector<double>& dat = data[idat];
    const double r2 = dat[0];
    if (r2 >= lower && r2 <= upper) {
      const double en = dat[1];
      double z = r2z_(r2, lower, upper);
      DEBUG("z " << MAX_PRECISION << z);
      if (z <= 0.) {
        ASSERT(z + NEAR_ZERO >= 0, "z:" << MAX_PRECISION << z);
        z = NEAR_ZERO;
      }
      if (z >= 1.) {
        ASSERT(z - NEAR_ZERO <= 1, "z:" << MAX_PRECISION << z);
        z = 1. - NEAR_ZERO;
      }
      DEBUG("z " << MAX_PRECISION << z);
      const double entab = table.linear_interpolation(z);
      const double crit = r2*std::abs(std::exp(-beta*en) - std::exp(-beta*entab));
      if (!filename.empty()) {
        file << r2 << " " << z << " " << en << " " << entab << " " << en - entab << " " << crit << std::endl;
      }
      if (crit > *max_criteria) {
        max_z = z;
        *max_criteria = crit;
      }
      //INFO("r2: " << r2 << " z: " << z << " en: " << en << " tab: " << entab << " diff: " << en - entab);
    }
  }
  ASSERT(*max_criteria > 0, "max_criteria: " << *max_criteria);
  INFO("max_z:" << max_z << " max_criteria:" << *max_criteria);
  const int max_bin = table.value_to_lowest_bin(max_z);
  INFO("max_bin:" << max_bin);
  return max_bin;
}

void RecursiveTablePotential::precompute(const Configuration& config) {
  Model::precompute(config);
  const ModelParams& existing = config.model_params();
  //const ModelParam& cutoff = existing.select("cutoff");
  t2index_.assign(existing.size(), -1);
  site_types_.clear();
  for (const std::string& sname : site_type_names_) {
    site_types_.push_back(config.site_type_name_to_index(sname));
  }
  for (int t1 = 0; t1 < static_cast<int>(site_types_.size()); ++t1) {
    const int type1 = site_types_[t1];
    t2index_[type1] = t1;
  }
  resize(1, 1, &cutoff_sq_);
  resize(1, 1, &inner_sq_);
  WARN("Assumes one site type.");
  t2index_.resize(1);
  t2index_[0] = 0;
  cutoff_sq_[0][0] = std::pow(config.model_params().select("cutoff").value(0), 2);
  if (energy_table_.size() == 0) {
    ASSERT(!mayer_training_file_.empty(), "requires mayer_training_file");
    std::ifstream file(mayer_training_file_);
    std::vector<std::vector<double> > data;
    std::string line;
    while (!file.eof()) {
      std::getline(file, line);
      //INFO("line: " << line);
      if (!line.empty()) {
        std::vector<std::string> dat = split(line, ',');
        std::vector<double> datd;
        //for (int sti = 0; sti < static_cast<int>(dat.size()); ++sti) {
        for (const std::string& st : dat) {
        //  const std::string& st = dat[sti];
          const double val = str_to_double(st);
          //if (sti == 0 && (val < inner_sq_[0][0] || val > cutoff_sq_[0][0])) {
          //  break;
          //}
          datd.push_back(val);
        }
        //INFO("datd:" << feasst_str(datd));
        if (datd.size() > 0) {
          data.push_back(datd);
        }
      }
    }
    DEBUG("sorting");
    std::sort(data.begin(), data.end());
    DEBUG("sorted");
    DEBUG("remove duplicants / redundant");
    auto it = std::unique(data.begin(), data.end());
    data.erase(it, data.end());

    // find the lowest distance with less than hard_limit in energy
    // remove data below that distance
    int found = -1;
    for (int index = 0; index < static_cast<int>(data.size()); ++index) {
      const std::vector<double>& en = data[index];
    //for (const std::vector<double>& en : data) {
      if (en[1] < hard_limit_u_) {
        inner_sq_[0][0] = en[0];
        found = index;
        DEBUG("found inner distance of " << en[0]);
        break;
      }
    }
    ASSERT(found != -1, "lowest inner distance not found");
    data.erase(data.begin(), data.begin() + found);

//    // print to screen
//    DEBUG(data.size());
//    for (const auto& d : data) {
//      DEBUG(feasst_str(d));
//    }

    resize(1, 1, &energy_table_);
    energy_table_[0][0] = std::make_unique<RecursiveTable1D>(argtype({{"num", str(size_)}}));
    RecursiveTable1D * tab = energy_table_[0][0].get();
    // find "size" equally spaced segments linearly interpolated from closest data
    std::vector<double> rsqs;
    build_table_(inner_sq_[0][0], cutoff_sq_[0][0], data, tab);
    int iteration = 0;
    double criteria;
    int max_bin = analyze_table_(inner_sq_[0][0], cutoff_sq_[0][0], data, *tab, beta_, "test"+str(iteration), &criteria);
    RecursiveTable1D nested(argtype({{"num", str(size_)}}));
    std::vector<int> bins = {max_bin};
    while (criteria > min_criteria_) {
      ++iteration;
      build_table_(z2r_(tab->bin_to_value(max_bin    ), inner_sq_[0][0], cutoff_sq_[0][0]),
                   z2r_(tab->bin_to_value(max_bin + 1), inner_sq_[0][0], cutoff_sq_[0][0]),
                   data, &nested);
      tab->insert(max_bin, nested);
      max_bin = analyze_table_(inner_sq_[0][0], cutoff_sq_[0][0], data, *tab, beta_, "test"+str(iteration), &criteria);
      ASSERT(iteration < 1e4, "Cannot get criteria:" << criteria <<
        " below min_criteria:" << min_criteria_ << " within iteration:"
        << iteration);
      if (criteria > min_criteria_) {
        ASSERT(!find_in_list(max_bin, bins), "At iteration:" << iteration <<
          ", max_bin:" << max_bin << " was already nested but found to have "
          << "criteria:" << criteria << ", which is higher than min_criteria:"
          << min_criteria_);
      }
      bins.push_back(max_bin);
    }
  }
}

double RecursiveTablePotential:: energy(
    const double squared_distance,
    const int type1,
    const int type2,
    const ModelParams& model_params) {
  TRACE("*** RecursiveTablePotential ***");
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
  const double inner_sq = inner_sq_[tabtype1][tabtype2];
  TRACE("squared_distance " << squared_distance);
  if (squared_distance < inner_sq) {
    TRACE("hard overlap");
    return NEAR_INFINITY;
  } else {
    const double rhsq = inner_sq_[tabtype1][tabtype2];
    TRACE("rhsq " << rhsq);
    const double rcsq = cutoff_sq_[tabtype1][tabtype2];
    TRACE("rcsq " << rcsq);
    const double z = (squared_distance - rhsq)/(rcsq - rhsq);
    TRACE("z " << z);
    ASSERT(z >= 0 && z <= 1, "z: " << z);
    TRACE("tab size " << energy_table_.size());
    TRACE("tab size " << energy_table_[0].size());
    const double en = energy_table_[tabtype1][tabtype2]->linear_interpolation(z);
    TRACE("en " << en);
    return en;
  }
}

FEASST_MAPPER(RecursiveTablePotential,);

RecursiveTablePotential::RecursiveTablePotential(std::istream& istr) : ModelTwoBody(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version >= 6765 && version <= 6765, "unrecognized version: " << version);
  feasst_deserialize(&inner_sq_, istr);
  feasst_deserialize(&cutoff_sq_, istr);
  feasst_deserialize(&site_types_, istr);
  feasst_deserialize(&site_type_names_, istr);
  feasst_deserialize(&t2index_, istr);
  feasst_deserialize(&energy_table_, istr);
  feasst_deserialize(&mayer_training_file_, istr);
  feasst_deserialize(&hard_limit_u_, istr);
  feasst_deserialize(&size_, istr);
  feasst_deserialize(&beta_, istr);
  feasst_deserialize(&min_criteria_, istr);
}

void RecursiveTablePotential::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_model_(ostr);
  feasst_serialize_version(6765, ostr);
  feasst_serialize(inner_sq_, ostr);
  feasst_serialize(cutoff_sq_, ostr);
  feasst_serialize(site_types_, ostr);
  feasst_serialize(site_type_names_, ostr);
  feasst_serialize(t2index_, ostr);
  feasst_serialize(energy_table_, ostr);
  feasst_serialize(mayer_training_file_, ostr);
  feasst_serialize(hard_limit_u_, ostr);
  feasst_serialize(size_, ostr);
  feasst_serialize(beta_, ostr);
  feasst_serialize(min_criteria_, ostr);
}

}  // namespace feasst
