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

std::shared_ptr<Model> Model::create(std::istream& istr) const {
  FATAL("not implemented");
}
std::shared_ptr<Model> Model::create(argtype * args) const {
  FATAL("not implemented");
}

Model::Model(std::istream& istr) {
  istr >> class_name_;
}

}  // namespace feasst
