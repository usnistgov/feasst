#include "utils/include/serialize.h"
#include "math/include/constants.h"  // PI
#include "math/include/random.h"
#include "math/include/utils_math.h"  // round
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
  const int angle_type = feasst::round(select->property("angle_type"));
  const Angle& angle = system->configuration().unique_type(
    select->particle_type()).angle(angle_type);
  angle_ = angle.property("theta0");
  if (angle.has_property("spring_constant")) {
    spring_constant_ = angle.property("spring_constant");
  }
}

double PerturbDistanceAngle::random_angle(Random * random,
    const double beta,
    const int dimension) const {
  if (std::abs(spring_constant_ + 1) < NEAR_ZERO) {
    return angle_;
  }
  double spring = spring_constant_/beta;
  if (std::abs(spring_constant_) < NEAR_ZERO) spring = 0.;
  return random->bond_angle(angle_, spring, 2, dimension);
}

void PerturbDistanceAngle::move(System * system,
                                TrialSelect * select,
                                Random * random) {
  if (origin_.dimension() == 0) {
    origin_.set_to_origin(system->configuration().dimension());
  }
  Select * mobile = select->get_mobile();
  Position * site = mobile->get_site_position(0, 0);
  DEBUG("mobile " << mobile->str());
  DEBUG("old pos " << site->str());

  /*
    For given sites j, k (anchors), place site, i, according to bond angle and
    length

    k - j
         \
   angle  i

    r_jk = r_j - r_k points from k to j
    - set the bond length by normalizing r_jk and multiplying by bond length
    - set the bond angle by rotating r_jk by 180-angle about axis orthogonal
      to r_jk.
   */
  // set site to the vector r_jk = r_j - r_k and store this vector
  const Position& rj = select->anchor_position(0, 0, *system);
  const Position& rk = select->anchor_position(0, 1, *system);
  rjk_ = rj;
  rjk_.subtract(rk);
  *site = rjk_;

  // normalize site position, centered at rj, and multiply by bond length.
  site->normalize();
  site->multiply(random_distance(random,
    system->thermo_params().beta(),
    system->dimension()));
  DEBUG("rjk " << rjk_.str());

  // rotate site by (PI-bond_angle) about vector orthogonal to r_jk
  orthogonal_jk_.orthogonal(*site);
  DEBUG("ortho " << orthogonal_jk_.str());
  const double ran_angle = random_angle(random, system->thermo_params().beta(),
    system->dimension());
  rot_mat_.axis_angle(orthogonal_jk_, PI - ran_angle);
  DEBUG("site == rj: " << site->str());
  rot_mat_.rotate(origin_, site);
  DEBUG("site rotated to angle: " << site->str());

  // randomly spin site about rjk.
  rot_mat_.axis_angle(rjk_, 2*PI*random->uniform());
  rot_mat_.rotate(origin_, site);

  site->add(rj);  // return to frame of reference

  DEBUG("new pos " << site->str());
  system->get_configuration()->update_positions(select->mobile());
}

PerturbDistanceAngle::PerturbDistanceAngle(std::istream& istr)
  : PerturbDistance(istr) {
  ASSERT(class_name_ == "PerturbDistanceAngle", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(788 == version, "mismatch version: " << version);
  feasst_deserialize(&angle_, istr);
  feasst_deserialize(&spring_constant_, istr);
}

void PerturbDistanceAngle::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_perturb_distance_(ostr);
  feasst_serialize_version(788, ostr);
  feasst_serialize(angle_, ostr);
  feasst_serialize(spring_constant_, ostr);
}

}  // namespace feasst
