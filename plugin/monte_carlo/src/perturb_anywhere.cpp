#include <cmath>
#include "utils/include/serialize.h"
#include "math/include/constants.h"
#include "configuration/include/domain.h"
#include "monte_carlo/include/perturb_anywhere.h"

namespace feasst {

PerturbAnywhere::PerturbAnywhere() {
  class_name_ = "PerturbAnywhere";
  rotate_.set_tunable(180.);
  rotate_.disable_tunable_();
  translate_.disable_tunable_();
  disable_tunable_();
}

class MapPerturbAnywhere {
 public:
  MapPerturbAnywhere() {
    auto obj = MakePerturbAnywhere();
    obj->deserialize_map()["PerturbAnywhere"] = obj;
  }
};

static MapPerturbAnywhere mapper_ = MapPerturbAnywhere();

void PerturbAnywhere::set_position(const Position& center,
                                   System * system,
                                   TrialSelect * select) {
  Position add = center;
  add.subtract(select->mobile().site_positions()[0][0]);
  translate_.move(add, system, select);
}

void PerturbAnywhere::move(System * system,
                           TrialSelect * select,
                           Random * random) {
  ASSERT(std::abs(rotate_.tunable().value() - 180.) < NEAR_ZERO,
    "rotation tunable should be 180");
  rotate_.move(system, select, random);
  system->configuration().domain().random_position(&random_in_box_, random);
  set_position(random_in_box_, system, select);
  DEBUG("anywhere: " << random_in_box_.str());
}

std::shared_ptr<Perturb> PerturbAnywhere::create(std::istream& istr) const {
  return std::make_shared<PerturbAnywhere>(istr);
}

PerturbAnywhere::PerturbAnywhere(std::istream& istr)
  : PerturbMove(istr) {
  ASSERT(class_name_ == "PerturbAnywhere", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(461 == version, "mismatch version: " << version);
  feasst_deserialize_fstobj(&rotate_, istr);
}

void PerturbAnywhere::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_perturb_(ostr);
  feasst_serialize_version(461, ostr);
  feasst_serialize_fstobj(rotate_, ostr);
}

}  // namespace feasst
