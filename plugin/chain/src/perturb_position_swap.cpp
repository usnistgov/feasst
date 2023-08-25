#include <cmath>
#include "utils/include/serialize.h"
#include "math/include/constants.h"
#include "math/include/utils_math.h"
#include "math/include/random.h"
#include "chain/include/perturb_position_swap.h"

namespace feasst {

PerturbPositionSwap::PerturbPositionSwap(argtype args) : PerturbPositionSwap(&args) {
  FEASST_CHECK_ALL_USED(args);
}
PerturbPositionSwap::PerturbPositionSwap(argtype * args) : PerturbMove(args) {
  class_name_ = "PerturbPositionSwap";
  disable_tunable_();
}

class MapPerturbPositionSwap {
 public:
  MapPerturbPositionSwap() {
    auto obj = MakePerturbPositionSwap();
    obj->deserialize_map()["PerturbPositionSwap"] = obj;
  }
};

static MapPerturbPositionSwap mapper_ = MapPerturbPositionSwap();

std::shared_ptr<Perturb> PerturbPositionSwap::create(std::istream& istr) const {
  return std::make_shared<PerturbPositionSwap>(istr);
}

void PerturbPositionSwap::move(const bool is_position_held,
                           System * system,
                           TrialSelect * select,
                           Random * random) {
  if (is_position_held) return;
  Select * mobile = select->get_mobile();
  Position * site0 = mobile->get_site_position(0, 0);
  Position * site1 = mobile->get_site_position(0, 1);
  INFO("mobile " << mobile->str());
  INFO("old pos " << site0->str() << " " << site1->str());
  Position site_tmp = *site0;
  *site0 = *site1;
  *site1 = site_tmp;
  INFO("new pos " << site0->str() << " " << site1->str());
  system->get_configuration()->update_positions(select->mobile());
}

PerturbPositionSwap::PerturbPositionSwap(std::istream& istr)
  : PerturbMove(istr) {
  // HWH can't check this if this is a base class
  // ASSERT(class_name_ == "PerturbPositionSwap", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(7841 == version, "mismatch version: " << version);
}

void PerturbPositionSwap::serialize_perturb_position_swap_(std::ostream& ostr) const {
  serialize_perturb_(ostr);
  feasst_serialize_version(7841, ostr);
}

void PerturbPositionSwap::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_perturb_position_swap_(ostr);
}

}  // namespace feasst
