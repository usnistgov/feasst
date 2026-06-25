#include "utils/include/serialize_extra.h"
#include "utils/include/arguments.h"
#include "math/include/position.h"
#include "configuration/include/configuration.h"
#include "configuration/include/model_params.h"
#include "system/include/model.h"

namespace feasst {

Model::Model() { class_name_ = "Model"; }
Model::~Model() {}

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

double Model::energy3body(
  const Position& r12,
  const Position& r13,
  const double squared_distance12,
  const double squared_distance13,
  const int type1,
  const int type2,
  const int type3,
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

void Model::precompute(Configuration * config, ModelParams * params) {
  epsilon_index_ = params->index("epsilon");
  sigma_index_ = params->index("sigma");
  cutoff_index_ = params->index("cutoff");
  charge_index_ = params->index("charge");
}

void Model::precompute(Configuration * config) {
  precompute(config, config->get_model_params());
}

}  // namespace feasst
