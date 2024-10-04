#include "system/include/dont_visit_model.h"
#include "utils/include/serialize.h"

namespace feasst {

FEASST_MAPPER(DontVisitModel,);

std::shared_ptr<VisitModel> DontVisitModel::create(std::istream& istr) const {
  return std::make_shared<DontVisitModel>(istr);
}

DontVisitModel::DontVisitModel(std::istream& istr) : VisitModel(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(104 == version, version);
}

void DontVisitModel::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_visit_model_(ostr);
  feasst_serialize_version(104, ostr);
}

}  // namespace feasst
