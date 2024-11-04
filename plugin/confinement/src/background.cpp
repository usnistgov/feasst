#include <cmath>
#include <vector>
#include "utils/include/arguments.h"
#include "utils/include/serialize.h"
#include "confinement/include/background.h"

namespace feasst {

Background::Background(argtype * args) : VisitModel(args) {
  class_name_ = "Background";
  constant_ = dble("constant", args);
}
Background::Background(argtype args) : Background(&args) {
  feasst_check_all_used(args);
}
Background::~Background() {}

void Background::compute(
    ModelOneBody * model,
    const ModelParams& model_params,
    const Select& selection,
    Configuration * config,
    const int group_index) {
  zero_energy();
}

FEASST_MAPPER(Background, argtype({{"constant", "0"}}));

std::shared_ptr<VisitModel> Background::create(std::istream& istr) const {
  return std::make_shared<Background>(istr);
}

std::shared_ptr<VisitModel> Background::create(argtype * args) const {
  return std::make_shared<Background>(args);
}

Background::Background(std::istream& istr) : VisitModel(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(3905 == version, version);
  feasst_deserialize(&constant_, istr);
}

void Background::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_visit_model_(ostr);
  feasst_serialize_version(3905, ostr);
  feasst_serialize(constant_, ostr);
}

void Background::compute(
    ModelOneBody * model,
    const ModelParams& model_params,
    Configuration * config,
    const int group_index) {
  set_energy(constant_);
}

}  // namespace feasst
