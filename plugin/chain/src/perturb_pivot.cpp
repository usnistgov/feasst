#include "utils/include/serialize.h"
#include "chain/include/perturb_pivot.h"
#include "math/include/random.h"

namespace feasst {

class MapPerturbPivot {
 public:
  MapPerturbPivot() {
    auto obj = MakePerturbPivot();
    obj->deserialize_map()["PerturbPivot"] = obj;
  }
};

static MapPerturbPivot mapper_ = MapPerturbPivot();

std::shared_ptr<Perturb> PerturbPivot::create(std::istream& istr) const {
  return std::make_shared<PerturbPivot>(istr);
}

PerturbPivot::PerturbPivot(std::istream& istr)
  : PerturbRotate(istr) {
  ASSERT(class_name_ == "PerturbPivot", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(445 == version, "mismatch version: " << version);
}

void PerturbPivot::serialize_perturb_pivot_(std::ostream& ostr) const {
  serialize_perturb_rotate_(ostr);
  feasst_serialize_version(445, ostr);
}

void PerturbPivot::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_perturb_pivot_(ostr);
}

void PerturbPivot::move(const bool is_position_held,
                        System * system,
                        TrialSelect * select,
                        Random * random) {
  if (is_position_held) return;
  const Position& pivot = select->anchor_position(0, 0, *system);
  DEBUG("piv " << pivot.str());
  PerturbRotate::move(system, select, random, pivot);
  DEBUG(select->mobile().site_positions()[0][0].str());
}
}  // namespace feasst
