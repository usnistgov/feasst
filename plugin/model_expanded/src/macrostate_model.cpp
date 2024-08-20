#include "utils/include/serialize.h"
#include "utils/include/arguments.h"
#include "math/include/histogram.h"
#include "model_expanded/include/macrostate_model.h"

namespace feasst {

MacrostateModel::MacrostateModel(const Histogram& histogram,
    argtype * args) : Macrostate(histogram, args) {
  class_name_ = "MacrostateModel";
  constraint_ = ConstrainModelIndex(args);
}
MacrostateModel::MacrostateModel(const Histogram& histogram,
    argtype args) : MacrostateModel(histogram, &args) {
  feasst_check_all_used(args);
}
MacrostateModel::MacrostateModel(argtype args) :
    MacrostateModel(Histogram(&args), &args) {
  feasst_check_all_used(args);
}

class MapMacrostateModel {
 public:
  MapMacrostateModel() {
    auto hist = MakeHistogram({{"width", "1"}, {"max", "1"}});
    MacrostateModel(*hist).deserialize_map()["MacrostateModel"] =
      MakeMacrostateModel(*hist);
  }
};

static MapMacrostateModel mapper_ = MapMacrostateModel();

std::shared_ptr<Macrostate> MacrostateModel::create(std::istream& istr) const {
  return std::make_shared<MacrostateModel>(istr);
}

std::shared_ptr<Macrostate> MacrostateModel::create(argtype * args) const {
  return std::make_shared<MacrostateModel>(args);
}

MacrostateModel::MacrostateModel(std::istream& istr)
  : Macrostate(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 9753, "version mismatch: " << version);
  feasst_deserialize_fstobj(&constraint_, istr);
}

void MacrostateModel::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_macrostate_(ostr);
  feasst_serialize_version(9753, ostr);
  feasst_serialize_fstobj(constraint_, ostr);
}

double MacrostateModel::value(const System& system,
    const Criteria& criteria,
    const Acceptance& acceptance) const {
  return constraint_.model_index(system, acceptance);
}

}  // namespace feasst
