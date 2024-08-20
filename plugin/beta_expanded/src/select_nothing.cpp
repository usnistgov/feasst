#include "utils/include/serialize.h"
#include "utils/include/arguments.h"
#include "math/include/random.h"
#include "beta_expanded/include/select_nothing.h"

namespace feasst {

SelectNothing::SelectNothing(argtype args) : SelectNothing(&args) {
  feasst_check_all_used(args);
}
SelectNothing::SelectNothing(argtype * args) : TrialSelect(args) {
  class_name_ = "SelectNothing";
}

class MapSelectNothing {
 public:
  MapSelectNothing() {
    auto obj = MakeSelectNothing();
    obj->deserialize_map()["SelectNothing"] = obj;
  }
};

static MapSelectNothing mapper_ = MapSelectNothing();

std::shared_ptr<TrialSelect> SelectNothing::create(std::istream& istr) const {
  return std::make_shared<SelectNothing>(istr);
}

SelectNothing::SelectNothing(std::istream& istr)
  : TrialSelect(istr) {
  // ASSERT(class_name_ == "SelectNothing", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(6947 == version, "mismatch version: " << version);
}

void SelectNothing::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_select_(ostr);
  feasst_serialize_version(6947, ostr);
}

}  // namespace feasst
