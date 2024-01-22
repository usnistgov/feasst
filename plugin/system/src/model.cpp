#include "configuration/include/model_params.h"
#include "system/include/model.h"
#include "utils/include/serialize.h"

namespace feasst {

std::map<std::string, std::shared_ptr<Model> >& Model::deserialize_map() {
  static std::map<std::string, std::shared_ptr<Model> >* ans =
     new std::map<std::string, std::shared_ptr<Model> >();
  return *ans;
}

std::shared_ptr<Model> Model::deserialize(std::istream& istr) {
  return template_deserialize(deserialize_map(), istr,
    // true argument denotes rewinding to reread class name
    // this allows derived class constructor to read class name.
    true);
}

std::shared_ptr<Model> Model::factory(const std::string name, argtype * args) {
  return template_factory(deserialize_map(), name, args);
}

double Model::energy(
  const double squared_distance,
  const int type1,
  const int type2,
  const ModelParams& model_params) { FATAL("not implemented"); }

void Model::serialize(std::ostream& ostr) const { FATAL("not implemented"); }

void Model::serialize_model_(std::ostream& ostr) const {
  feasst_serialize_version(2094, ostr);
  feasst_serialize(epsilon_index_, ostr);
  feasst_serialize(sigma_index_, ostr);
  feasst_serialize(cutoff_index_, ostr);
  feasst_serialize(charge_index_, ostr);
}

std::shared_ptr<Model> Model::create(std::istream& istr) const {
  FATAL("not implemented");
}
std::shared_ptr<Model> Model::create(argtype * args) const {
  FATAL("not implemented");
}

int Model::model_index() const {
  FATAL("not implemented");
}

void Model::set_model_index(const int index) {
  FATAL("not implemented");
}

Model::Model(std::istream& istr) {
  istr >> class_name_;
  const int version = feasst_deserialize_version(istr);
  ASSERT(2094 == version, "version mismatch: " << version);
  feasst_deserialize(&epsilon_index_, istr);
  feasst_deserialize(&sigma_index_, istr);
  feasst_deserialize(&cutoff_index_, istr);
  feasst_deserialize(&charge_index_, istr);
}

void Model::precompute(const ModelParams& existing) {
  epsilon_index_ = existing.index("epsilon");
  sigma_index_ = existing.index("sigma");
  cutoff_index_ = existing.index("cutoff");
  charge_index_ = existing.index("charge");
}

}  // namespace feasst
