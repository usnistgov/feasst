#include "utils/include/serialize.h"
#include "math/include/constants.h"  // PI
#include "math/include/random.h"
#include "monte_carlo/include/perturb_distance_angle.h"

namespace feasst {

PerturbDistanceAngle::PerturbDistanceAngle(const argtype& args)
  : PerturbDistance(args) {
  class_name_ = "PerturbDistanceAngle";
}

class MapPerturbDistanceAngle {
 public:
  MapPerturbDistanceAngle() {
    auto obj = MakePerturbDistanceAngle();
    obj->deserialize_map()["PerturbDistanceAngle"] = obj;
  }
};

static MapPerturbDistanceAngle mapper_ = MapPerturbDistanceAngle();

std::shared_ptr<Perturb> PerturbDistanceAngle::create(std::istream& istr) const {
  return std::make_shared<PerturbDistanceAngle>(istr);
}

void PerturbDistanceAngle::precompute(TrialSelect * select, System * system) {
  PerturbDistance::precompute(select, system);
  angle_ = select->property("theta0");
  origin_.set_to_origin(system->configuration().dimension());
  rjk_ = origin_;
}

void PerturbDistanceAngle::move(System * system,
                                TrialSelect * select,
                                Random * random) {
  Select * mobile = select->get_mobile();
  Position * site = mobile->get_site_position(0, 0);
  DEBUG("mobile " << mobile->str());
  DEBUG("old pos " << site->str());

  // set site to the vector |r_j - r_k| and store this unit vector
  const Position& rj = select->anchor_position(0, 0, *system);
  const Position& rk = select->anchor_position(0, 1, *system);
  *site = rj;
  site->subtract(rk);
  rjk_ = *site;
  DEBUG("rjk " << rjk_.str());
  //rjk_.normalize();

  // rotate site by (PI-theta0) about vector orthogonal to r_jk
  orthogonal_jk_.orthogonal(*site);
  DEBUG("ortho " << orthogonal_jk_.str());
  rot_mat_.axis_angle(orthogonal_jk_, PI - angle_);
  DEBUG("site == rj: " << site->str());
  rot_mat_.rotate(origin_, site);
  DEBUG("site rotated to angle: " << site->str());

  // randomly spin site about rjk.
  rot_mat_.axis_angle(rjk_, 2*PI*random->uniform());
  rot_mat_.rotate(origin_, site);

  site->add(rj);  // return frame of reference

  DEBUG("new pos " << site->str());
  system->get_configuration()->update_positions(select->mobile());
}

PerturbDistanceAngle::PerturbDistanceAngle(std::istream& istr)
  : PerturbDistance(istr) {
  ASSERT(class_name_ == "PerturbDistanceAngle", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(788 == version, "mismatch version: " << version);
  feasst_deserialize(&angle_, istr);
}

void PerturbDistanceAngle::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_perturb_distance_(ostr);
  feasst_serialize_version(788, ostr);
  feasst_serialize(angle_, ostr);
}

}  // namespace feasst
