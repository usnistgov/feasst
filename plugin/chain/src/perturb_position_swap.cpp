#include <cmath>
#include "utils/include/arguments.h"
#include "utils/include/serialize.h"
#include "math/include/constants.h"
#include "math/include/utils_math.h"
#include "math/include/random.h"
#include "configuration/include/select.h"
#include "configuration/include/particle_factory.h"
#include "configuration/include/configuration.h"
#include "monte_carlo/include/trial_select.h"
#include "chain/include/perturb_position_swap.h"

namespace feasst {

PerturbPositionSwap::PerturbPositionSwap(argtype args) : PerturbPositionSwap(&args) {
  feasst_check_all_used(args);
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
                           Random * random,
                           Acceptance * acceptance) {
  if (is_position_held) return;
  Select * mobile = select->get_mobile();
  Position * site0 = mobile->get_site_position(0, 0);
  Position * site1;
  if (mobile->num_particles() == 1) {
    site1 = mobile->get_site_position(0, 1);
  } else if (mobile->num_particles() == 2) {
    site1 = mobile->get_site_position(1, 0);
  } else {
    FATAL("unrecognized number of mobile parts: " << mobile->num_particles());
  }
  DEBUG("mobile " << mobile->str());
  DEBUG("old pos " << site0->str() << " " << site1->str());
  Position site_tmp = *site0;
  *site0 = *site1;
  *site1 = site_tmp;
  DEBUG("new pos " << site0->str() << " " << site1->str());
  select->get_configuration(system)->update_positions(select->mobile());
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
