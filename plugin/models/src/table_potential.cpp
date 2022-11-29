#include <cmath>
#include <string>
#include <fstream>
#include "utils/include/serialize.h"
#include "math/include/utils_math.h"
#include "configuration/include/domain.h"
#include "configuration/include/configuration.h"
#include "system/include/model_two_body.h"
#include "models/include/table_potential.h"

namespace feasst {

TablePotential::TablePotential(argtype * args) : VisitModelInner(args) {
  class_name_ = "TablePotential";
  const std::string table_file = str("table_file", args, "-1");
  if (table_file != "-1") {
    read_table_(table_file);
  }
}
TablePotential::TablePotential(argtype args) : TablePotential(&args) {
  FEASST_CHECK_ALL_USED(args);
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
  file >> descript >> int_val;
  ASSERT(descript == "site_types", "format error: " << descript);
  const int num_sites = int_val;
  DEBUG("num_sites " << num_sites);

  // size arrays
  site_types_.resize(num_sites);
  resize(num_sites, num_sites, &inner_);
  resize(num_sites, num_sites, &inner_g_);
  resize(num_sites, num_sites, &cutoff_g_);
  resize(num_sites, num_sites, &energy_table_);
//  resize(num_sites, num_sites, &gamma_);
  for (int type = 0; type < num_sites; ++type) {
    file >> int_val;
    DEBUG("site " << int_val);
    site_types_[type] = int_val;
  }

  for (int itype = 0; itype < num_sites; ++itype) {
    for (int jtype = itype; jtype < num_sites; ++jtype) {
//      // read gamma
//      file >> descript >> double_val;
//      ASSERT(descript == "gamma", "format error: " << descript);
//      const double gamma = double_val;
//      DEBUG("gamma " << gamma);
//      gamma_[itype][jtype] = gamma;

      // read inner
      file >> descript >> double_val;
      ASSERT(descript == "inner", "format error: " << descript);
      const double inner = double_val;
      DEBUG("inner " << inner);
      ASSERT(inner > 0, "inner: " << inner << " must be > 0");
      inner_[itype][jtype] = inner;
      inner_g_[itype][jtype] = std::pow(inner, gamma_);

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

void TablePotential::precompute(Configuration * config) {
  VisitModelInner::precompute(config);
  const ModelParam& cutoff = config->model_params().select("cutoff");
  for (int t1 = 0; t1 < static_cast<int>(site_types_.size()); ++t1) {
    const int type1 = site_types_[t1];
    for (int t2 = 0; t2 < static_cast<int>(site_types_.size()); ++t2) {
      const int type2 = site_types_[t2];
      const double rc = cutoff.mixed_value(type1, type2);
      cutoff_g_[t1][t2] = std::pow(rc, gamma_);
    }
  }
}

void TablePotential::compute(
    const int part1_index,
    const int site1_index,
    const int part2_index,
    const int site2_index,
    const Configuration * config,
    const ModelParams& model_params,
    ModelTwoBody * model,
    const bool is_old_config,
    Position * relative,
    Position * pbc) {
  DEBUG("*** TablePotential ***");
  const Particle& part1 = config->select_particle(part1_index);
  const Site& site1 = part1.site(site1_index);
  const Particle& part2 = config->select_particle(part2_index);
  const Site& site2 = part2.site(site2_index);
  clear_ixn(part1_index, site1_index, part2_index, site2_index);
  int type1 = site1.type();
  int type2 = site2.type();

  // check if sites are within the global cutoff
  const double cutoff = model_params.select(cutoff_index()).mixed_values()[type1][type2];
  double squared_distance;
  config->domain().wrap_opt(site1.position(), site2.position(), relative, pbc, &squared_distance);
  DEBUG("squared_distance " << squared_distance);
  DEBUG("relative " << relative->str());
  if (squared_distance > cutoff*cutoff) {
    return;
  }
  DEBUG("inside global cut");

  // enforce type1 <= type2 to avoid redundant tables.
  bool flip = false;
  if (type1 > type2) {
    flip = true;
  }
  if (flip) {
    swap(&type1, &type2);
  }
  DEBUG("flip " << flip);

  // check the inner cutoff.
  const double inner = inner_[type1][type2];
  double en;
  if (squared_distance < inner*inner) {
    en = NEAR_INFINITY;
    DEBUG("hard overlap");
  } else {
    // const double gamma = gamma_[type1][type2];
    // DEBUG("gamma " << gamma);
    const double rhg = inner_g_[type1][type2];
    DEBUG("rhg " << rhg);
    const double rcg = cutoff_g_[type1][type2];
    DEBUG("rcg " << rcg);
    const double rg = 1./squared_distance;
    DEBUG("rg " << rg);
    const double z = (rg - rhg)/(rcg - rhg);
    DEBUG("z " << z);
    ASSERT(z >= 0 && z <= 1, "z: " << z);
    DEBUG("type1 " << type1 << " type2 " << type2);
    DEBUG("tab size " << energy_table_.size());
    DEBUG("tab size " << energy_table_[0].size());
    en = energy_table_[type1][type2].forward_difference_interpolation(z);
  }
  DEBUG("en " << en);
  update_ixn(en, part1_index, site1_index, type1, part2_index,
             site2_index, type2, squared_distance, pbc, is_old_config, *config);
}

class MapTablePotential {
 public:
  MapTablePotential() {
    auto obj = MakeTablePotential();
    obj->deserialize_map()["TablePotential"] = obj;
  }
};

static MapTablePotential mapper_ = MapTablePotential();

TablePotential::TablePotential(std::istream& istr) : VisitModelInner(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 8965, "unrecognized version: " << version);
  feasst_deserialize(&inner_, istr);
  feasst_deserialize(&inner_g_, istr);
  feasst_deserialize(&cutoff_g_, istr);
  feasst_deserialize(&site_types_, istr);
  feasst_deserialize_fstobj(&energy_table_, istr);
}

void TablePotential::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_visit_model_inner_(ostr);
  feasst_serialize_version(8965, ostr);
  feasst_serialize(inner_, ostr);
  feasst_serialize(inner_g_, ostr);
  feasst_serialize(cutoff_g_, ostr);
  feasst_serialize(site_types_, ostr);
  feasst_serialize_fstobj(energy_table_, ostr);
}

}  // namespace feasst
