#include "utils/include/serialize.h"
#include "math/include/constants.h"  // PI
#include "math/include/random.h"
#include "math/include/utils_math.h"  // round
#include "monte_carlo/include/perturb_distance_angle.h"

namespace feasst {

PerturbDistanceAngle::PerturbDistanceAngle(argtype args)
  : PerturbDistanceAngle(&args) {
  check_all_used(args);
}
PerturbDistanceAngle::PerturbDistanceAngle(argtype * args)
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
  ASSERT(select->has_property("angle_type"), "cannot obtain angle properties");
  angle_type_ = feasst::round(select->property("angle_type"));
  DEBUG("angle_type_ " << angle_type_);
}

double PerturbDistanceAngle::random_angle_radians(const System& system,
    const TrialSelect* select,
    Random * random) {
  FATAL("not implemented");
}
double PerturbDistanceAngle::random_angle_radians(const System& system,
    const TrialSelect* select, Random * random, double * bond_energy) {
  const Angle& angle = system.configuration().unique_type(
    select->particle_type()).angle(angle_type_);
  const double beta = system.thermo_params().beta();
  ASSERT(angle_.deserialize_map().count(angle.model()) == 1,
    angle.model() << " not found");
  const BondThreeBody * model = angle_.deserialize_map()[angle.model()].get();
  const double radians = model->random_angle_radians(angle, beta, system.dimension(), random);
  *bond_energy += model->energy(radians, angle);
  DEBUG("bond_energy " << *bond_energy);
  return radians;
}

void PerturbDistanceAngle::move(System * system,
    TrialSelect * select,
    Random * random) {
  DEBUG(class_name());
  double bond_energy = 0.;
  const double distance = random_distance(*system, select, random, &bond_energy);
  DEBUG("bond_energy " << bond_energy);
  const double angle = random_angle_radians(*system, select, random, &bond_energy);
  DEBUG("final angle pert bond_energy " << bond_energy);
  select->add_exclude_energy(bond_energy);
  DEBUG("angle: " << angle);
  place_in_circle(distance, angle, system, select, random);
}

void PerturbDistanceAngle::place_in_circle(const double distance,
  const double angle,
  System * system,
  TrialSelect * select,
  Random * random) {
  const int dimension = system->configuration().dimension();
  if (origin_.dimension() == 0) origin_.set_to_origin(dimension);
  if (orthogonal_jk_.dimension() == 0) orthogonal_jk_.set_to_origin(dimension);
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
  DEBUG("rj: " << rj.str());
  const Position& rk = select->anchor_position(0, 1, *system);
  DEBUG("rk: " << rk.str());
  rjk_ = rj;
  rjk_.subtract(rk);
  *site = rjk_;
  DEBUG("rjk: " << rjk_.str());

  // normalize site position, centered at rj, and multiply by bond length.
  site->normalize();
  site->multiply(distance);
  DEBUG("rjk norm*L: " << rjk_.str());

  // rotate site by (PI-bond_angle). If 3D, about vector orthogonal to r_jk.
  if (dimension == 3) {
    orthogonal_jk_.orthogonal(*site);
    DEBUG("ortho " << orthogonal_jk_.str());
    rot_mat_.axis_angle(orthogonal_jk_, radians_to_degrees(PI - angle));
  } else {
    rot_mat_.axis_angle(orthogonal_jk_, radians_to_degrees(PI + angle));
  }
  DEBUG("site == rjk: " << site->str());
  rot_mat_.rotate(origin_, site);
  DEBUG("site rotated to angle: " << site->str());

  // If 3D, randomly spin site about rjk.
  if (dimension == 3) {
    rot_mat_.axis_angle(rjk_, 360.*random->uniform());
    rot_mat_.rotate(origin_, site);
  }

  site->add(rj);  // return to frame of reference

  DEBUG("new pos " << site->str());
  system->get_configuration()->update_positions(select->mobile());
}

PerturbDistanceAngle::PerturbDistanceAngle(std::istream& istr)
  : PerturbDistance(istr) {
  //ASSERT(class_name_ == "PerturbDistanceAngle", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(788 == version, "mismatch version: " << version);
  feasst_deserialize(&angle_type_, istr);
}

void PerturbDistanceAngle::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_perturb_distance_angle_(ostr);
}

void PerturbDistanceAngle::serialize_perturb_distance_angle_(
    std::ostream& ostr) const {
  serialize_perturb_distance_(ostr);
  feasst_serialize_version(788, ostr);
  feasst_serialize(angle_type_, ostr);
}

}  // namespace feasst
