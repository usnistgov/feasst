#include "utils/include/serialize.h"
#include "monte_carlo/include/perturb_distance.h"
#include "math/include/random.h"

namespace feasst {

PerturbDistance::PerturbDistance(const argtype& args) : PerturbMove(args) {
  class_name_ = "PerturbDistance";
  disable_tunable_();
}

class MapPerturbDistance {
 public:
  MapPerturbDistance() {
    auto obj = MakePerturbDistance();
    obj->deserialize_map()["PerturbDistance"] = obj;
  }
};

static MapPerturbDistance mapper_ = MapPerturbDistance();

std::shared_ptr<Perturb> PerturbDistance::create(std::istream& istr) const {
  return std::make_shared<PerturbDistance>(istr);
}

void PerturbDistance::precompute(TrialSelect * select, System * system) {
  // determine the bond length
  // or input the bond length
  if (select->has_property("bond_length")) {
    distance_ = select->property("bond_length");
  } else {
    WARN("using default distance (typically for reptation): " << distance_);
  }
}

void PerturbDistance::move(System * system,
                           TrialSelect * select,
                           Random * random) {
  Select * mobile = select->get_mobile();
  Position * site = mobile->get_site_position(0, 0);
  DEBUG("mobile " << mobile->str());
  DEBUG("old pos " << site->str());
  random->unit_sphere_surface(site);
  site->multiply(distance_);
  site->add(select->anchor_position(0, 0, system));
  DEBUG("new pos " << site->str());
  system->get_configuration()->update_positions(select->mobile());
}

PerturbDistance::PerturbDistance(std::istream& istr)
  : PerturbMove(istr) {
  // HWH can't check this if this is a base class
  // ASSERT(class_name_ == "PerturbDistance", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(228 == version, "mismatch version: " << version);
  feasst_deserialize(&distance_, istr);
}

void PerturbDistance::serialize_perturb_distance_(std::ostream& ostr) const {
  serialize_perturb_(ostr);
  feasst_serialize_version(228, ostr);
  feasst_serialize(distance_, ostr);
}

void PerturbDistance::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_perturb_distance_(ostr);
}

}  // namespace feasst
