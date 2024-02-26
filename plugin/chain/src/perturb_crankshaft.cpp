#include "utils/include/debug.h"
#include "utils/include/serialize.h"
#include "chain/include/perturb_crankshaft.h"
#include "math/include/random.h"

namespace feasst {

class MapPerturbCrankshaft {
 public:
  MapPerturbCrankshaft() {
    auto obj = MakePerturbCrankshaft();
    obj->deserialize_map()["PerturbCrankshaft"] = obj;
  }
};

static MapPerturbCrankshaft mapper_ = MapPerturbCrankshaft();

std::shared_ptr<Perturb> PerturbCrankshaft::create(std::istream& istr) const {
  return std::make_shared<PerturbCrankshaft>(istr);
}

PerturbCrankshaft::PerturbCrankshaft(std::istream& istr)
  : PerturbRotate(istr) {
  ASSERT(class_name_ == "PerturbCrankshaft", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(555 == version, "mismatch version: " << version);
}

void PerturbCrankshaft::serialize_perturb_crankshaft_(std::ostream& ostr) const {
  serialize_perturb_rotate_(ostr);
  feasst_serialize_version(555, ostr);
}

void PerturbCrankshaft::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_perturb_crankshaft_(ostr);
}

void PerturbCrankshaft::move(const bool is_position_held,
    System * system,
    TrialSelect * select,
    Random * random,
    Acceptance * acceptance) {
  if (is_position_held) return;
  const Position& pivot = select->mobile().site_positions()[0].front();
  axis_ = select->mobile().site_positions()[0].back();
  axis_.subtract(pivot);
  axis_.normalize();
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
