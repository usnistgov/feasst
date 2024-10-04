#include <cmath>
#include "utils/include/io.h"
#include "utils/include/arguments.h"
#include "utils/include/utils.h"
#include "utils/include/serialize.h"
#include "math/include/constants.h"
#include "math/include/table.h"
#include "configuration/include/model_params.h"
#include "system/include/model_two_body_table.h"

namespace feasst {

FEASST_MAPPER(ModelTwoBodyTable,);

ModelTwoBodyTable::ModelTwoBodyTable(argtype args) : ModelTwoBodyTable(&args) {
  feasst_check_all_used(args);
}

ModelTwoBodyTable::ModelTwoBodyTable(argtype * args) {
  class_name_ = "ModelTwoBodyTable";
  const double hard_sphere_threshold = dble("hard_sphere_threshold", args, 0.85);
  hard_sphere_threshold_inv_sq_ = 1/std::pow(hard_sphere_threshold, 2);
  cutoff_inv_sq_ = std::make_unique<CutOff>();
}
ModelTwoBodyTable::~ModelTwoBodyTable() {}

void ModelTwoBodyTable::resize(const int num_site_types) {
  feasst::resize(num_site_types, num_site_types, &table_);
}

void ModelTwoBodyTable::precompute(const ModelParams& existing) {
  Model::precompute(existing);
  if (cutoff_inv_sq_->size() == 0) {
    const ModelParam& rc = existing.select("cutoff");
    for (int type1 = 0; type1 < existing.size(); ++type1) {
      cutoff_inv_sq_->add(std::pow(rc.value(type1), -2));
    }
    cutoff_inv_sq_->mix();
    for (int type1 = 0; type1 < existing.size(); ++type1) {
      for (int type2 = 0; type2 < existing.size(); ++type2) {
        cutoff_inv_sq_->set_mixed(type1, type2,
                                 std::pow(rc.mixed_value(type1, type2), -2));
      }
    }
  }
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
  INFO("finding sigma");
  const ModelParam& sig = model_params.select("sigma");
  INFO("num_types " << num_types);
  for (int type1 = 0; type1 < num_types; ++type1) {
    DEBUG("type1 " << type1);
    //for (int type2 = 0; type2 < num_types; ++type2) {
    for (int type2 = type1; type2 < num_types; ++type2) {
      DEBUG("type2 " << type2);
      const double sigma = sig.mixed_value(type1, type2);
      const double rhg = hard_sphere_threshold_inv_sq_/sigma/sigma;
      const double rcg = cutoff_inv_sq_->mixed_value(type1, type2);
      auto table = MakeTable1D({{"num", str(size)}});
      for (int bin = 0; bin < size; ++bin) {
        const double z = table->bin_to_value(bin);
        const double r = std::pow(z*(rcg - rhg) + rhg, -0.5);
        const double en = model->energy(r*r, type1, type2, model_params);
        DEBUG("z " << z << " r " << r << " en " << en);
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
  serialize_model_(ostr);
  feasst_serialize_version(5937, ostr);
  feasst_serialize(hard_sphere_threshold_inv_sq_, ostr);
  feasst_serialize(table_, ostr);
  feasst_serialize(cutoff_inv_sq_, ostr);
}

ModelTwoBodyTable::ModelTwoBodyTable(std::istream& istr) : ModelTwoBody(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(5937 == version, version);
  feasst_deserialize(&hard_sphere_threshold_inv_sq_, istr);
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
  feasst_deserialize(cutoff_inv_sq_, istr);
}

double ModelTwoBodyTable::energy(
    const double squared_distance,
    const int type1,
    const int type2,
    const ModelParams& model_params) {
  TRACE("squared_distance " << squared_distance);
  const double sigma =
    model_params.select(sigma_index()).mixed_value(type1, type2);
  const double rhg = hard_sphere_threshold_inv_sq_/sigma/sigma;
  TRACE("rhg " << rhg);
  TRACE("type1 " << type1);
  TRACE("type2 " << type2);
  TRACE("table size: " << table_.size());
  TRACE("table size 2: " << table_[type1].size());
  TRACE("table? " << table_[type1][type2]);
  if (1./squared_distance > rhg) {
    return NEAR_INFINITY;
  } else if (table_[type1][type2]) {
    const double rcg = cutoff_inv_sq_->mixed_value(type1, type2);
    TRACE("rcg " << rcg);
    const double z = (1./squared_distance - rhg)/(rcg - rhg);
    TRACE("z " << z);
    const double en = table_[type1][type2]->linear_interpolation(z);
    TRACE("en " << en);
    return en;
  } else {
    return 0.;
  }
}

}  // namespace feasst
