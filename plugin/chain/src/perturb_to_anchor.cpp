#include <cmath>
#include "utils/include/serialize.h"
#include "math/include/constants.h"
#include "math/include/utils_math.h"
#include "math/include/random.h"
#include "chain/include/perturb_to_anchor.h"

namespace feasst {

PerturbToAnchor::PerturbToAnchor(argtype args) : PerturbToAnchor(&args) {
  FEASST_CHECK_ALL_USED(args);
}
PerturbToAnchor::PerturbToAnchor(argtype * args) : PerturbMove(args) {
  class_name_ = "PerturbToAnchor";
  disable_tunable_();
}

class MapPerturbToAnchor {
 public:
  MapPerturbToAnchor() {
    auto obj = MakePerturbToAnchor();
    obj->deserialize_map()["PerturbToAnchor"] = obj;
  }
};

static MapPerturbToAnchor mapper_ = MapPerturbToAnchor();

std::shared_ptr<Perturb> PerturbToAnchor::create(std::istream& istr) const {
  return std::make_shared<PerturbToAnchor>(istr);
}

void PerturbToAnchor::move(const bool is_position_held,
                           System * system,
                           TrialSelect * select,
                           Random * random) {
  if (is_position_held) return;
  Select * mobile = select->get_mobile();
  Position * site = mobile->get_site_position(0, 0);
  DEBUG("mobile " << mobile->str());
  DEBUG("old pos " << site->str());
  *site = select->anchor_position(0, 0, *system);
  DEBUG("new pos " << site->str());
  select->get_configuration(system)->update_positions(select->mobile());
}

PerturbToAnchor::PerturbToAnchor(std::istream& istr)
  : PerturbMove(istr) {
  // HWH can't check this if this is a base class
  // ASSERT(class_name_ == "PerturbToAnchor", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(6506 == version, "mismatch version: " << version);
}

void PerturbToAnchor::serialize_perturb_to_anchor_(std::ostream& ostr) const {
  serialize_perturb_(ostr);
  feasst_serialize_version(6506, ostr);
}

void PerturbToAnchor::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_perturb_to_anchor_(ostr);
}

}  // namespace feasst
