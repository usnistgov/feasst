#include <cmath>
#include "utils/include/utils.h"
#include "utils/include/serialize.h"
#include "math/include/constants.h"
#include "math/include/table.h"
#include "system/include/model_two_body_table.h"

namespace feasst {

ModelTwoBodyTable::ModelTwoBodyTable(argtype args) : ModelTwoBodyTable(&args) {
  check_all_used(args);
}

ModelTwoBodyTable::ModelTwoBodyTable(argtype * args) {
  class_name_ = "ModelTwoBodyTable";
  hard_sphere_threshold_ = dble("hard_sphere_threshold", args, 0.2);
}

class MapModelTwoBodyTable {
 public:
  MapModelTwoBodyTable() {
    ModelTwoBodyTable().deserialize_map()["ModelTwoBodyTable"] =
      MakeModelTwoBodyTable();
  }
};

static MapModelTwoBodyTable map_model_two_body_table_ = MapModelTwoBodyTable();

void ModelTwoBodyTable::resize(const int num_site_types) {
  feasst::resize(num_site_types, num_site_types, &table_);
}

void ModelTwoBodyTable::set(std::shared_ptr<Table1D> table,
    const int type1,
    const int type2) {
  DEBUG("setting table " << type1 << " " << type2 << " " << table);
  table_[type1][type2] = table;
  if (type1 < type2) set(table, type2, type1);
//  // table exists?
//  for (int type1 = 0; type1 < static_cast<int>(table_.size()); ++type1) {
//    for (int type2 = 0; type2 < static_cast<int>(table_.size()); ++type2) {
//      INFO("tab " << type1 << " " << type2 << " " << table_[type1][type2]);
//    }
//  }
}

void ModelTwoBodyTable::set(const ModelParams& model_params,
    const int size,
    const int num_types,
    Model * model) {
  resize(num_types);
  const double rh = hard_sphere_threshold_;
  DEBUG("num_types " << num_types);
  for (int type1 = 0; type1 < num_types; ++type1) {
    //INFO("type1 " << type1);
    //for (int type2 = 0; type2 < num_types; ++type2) {
    for (int type2 = type1; type2 < num_types; ++type2) {
      //INFO("type2 " << type2);
      const double rc = model_params.mixed_cutoff()[type1][type2];
      auto table = MakeTable1D({{"num", str(size)}});
      for (int bin = 0; bin < size; ++bin) {
        const double z = table->bin_to_value(bin);
        const double r = z*(rc - rh) + rh;
        const double en = model->energy(r*r, type1, type2, model_params);
        //INFO("z " << z << " r " << r << " en " << en);
        table->set_data(bin, en);
        //table->set_data(bin, model->energy(r*r, type1, type2, model_params));
      }
      set(table, type1, type2);
    }
  }
}

void ModelTwoBodyTable::set(const ModelParams& model_params,
    const int size,
    const int num_types,
    std::shared_ptr<Model> model) {
  set(model_params, size, num_types, model.get());
}

void ModelTwoBodyTable::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_model_two_body_table_(ostr);
}

void ModelTwoBodyTable::serialize_model_two_body_table_(
    std::ostream& ostr) const {
  feasst_serialize_version(5937, ostr);
  feasst_serialize(hard_sphere_threshold_, ostr);
  feasst_serialize(table_, ostr);
}

ModelTwoBodyTable::ModelTwoBodyTable(std::istream& istr) : ModelTwoBody(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(5937 == version, version);
  feasst_deserialize(&hard_sphere_threshold_, istr);
  // HWH for unknown reasons, this function template does not work.
  //feasst_deserialize(&table_, istr);
  { int dim1;
    istr >> dim1;
    table_.resize(dim1);
    for (int index1 = 0; index1 < dim1; ++index1) {
      int dim2;
      istr >> dim2;
      table_[index1].resize(dim2);
      for (int index2 = 0; index2 < dim2; ++index2) {
        int existing;
        istr >> existing;
        if (existing != 0) {
          table_[index1][index2] = std::make_shared<Table1D>(istr);
        }
      }
    }
  }
}

double ModelTwoBodyTable::energy(
    const double squared_distance,
    const int type1,
    const int type2,
    const ModelParams& model_params) {
  const double distance = std::sqrt(squared_distance);
  TRACE("distance " << distance);
  const double rh = hard_sphere_threshold_;
  TRACE("rh " << rh);
  TRACE("type1 " << type1);
  TRACE("type2 " << type2);
  TRACE("table size: " << table_.size());
  TRACE("table size 2: " << table_[type1].size());
  TRACE("table? " << table_[type1][type2]);
  if (distance < rh) {
    return NEAR_INFINITY;
  } else if (table_[type1][type2]) {
    const double cutoff = model_params.mixed_cutoff()[type1][type2];
    TRACE("cutoff " << cutoff);
    const double z = (distance - rh)/(cutoff - rh);
    TRACE("z " << z);
    const double en = table_[type1][type2]->linear_interpolation(z);
    TRACE("en " << en);
    return en;
  } else {
    return 0.;
  }
}

}  // namespace feasst
