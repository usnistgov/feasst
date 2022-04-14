#include "utils/include/serialize.h"
#include "chain/include/perturb_crankshaft_small.h"
#include "math/include/random.h"

namespace feasst {

class MapPerturbCrankshaftSmall {
 public:
  MapPerturbCrankshaftSmall() {
    auto obj = MakePerturbCrankshaftSmall();
    obj->deserialize_map()["PerturbCrankshaftSmall"] = obj;
  }
};

static MapPerturbCrankshaftSmall mapper_ = MapPerturbCrankshaftSmall();

std::shared_ptr<Perturb> PerturbCrankshaftSmall::create(std::istream& istr) const {
  return std::make_shared<PerturbCrankshaftSmall>(istr);
}

PerturbCrankshaftSmall::PerturbCrankshaftSmall(std::istream& istr)
  : PerturbRotate(istr) {
  ASSERT(class_name_ == "PerturbCrankshaftSmall", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(555 == version, "mismatch version: " << version);
}

void PerturbCrankshaftSmall::serialize_perturb_crankshaft_(std::ostream& ostr) const {
  serialize_perturb_rotate_(ostr);
  feasst_serialize_version(555, ostr);
}

void PerturbCrankshaftSmall::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_perturb_crankshaft_(ostr);
}

void PerturbCrankshaftSmall::move(const bool is_position_held,
    System * system, TrialSelect * select,
    Random * random) {
  if (is_position_held) return;
  DEBUG("anchor " << select->anchor().str());
  const Position& pivot = select->anchor_position(0, 0, *system);
  DEBUG("anchor0 " << pivot.str())
  axis_ = select->anchor_position(0, 1, *system);
  DEBUG("anchor1 " << axis_.str())
  axis_.subtract(pivot);
  axis_.normalize();
  DEBUG("axis " << axis_.str())
  const double max_angle = tunable().value();
  const double angle = random->uniform_real(-max_angle, max_angle);
  if (rot_mat_.num_rows() == 0) {
    rot_mat_.set_size(axis_.size(), axis_.size());
  }
  rot_mat_.axis_angle_opt(axis_, angle);
  DEBUG("mobile " << select->mobile().str());
  PerturbRotate::move(pivot, rot_mat_, system, select);
}

}  // namespace feasst
