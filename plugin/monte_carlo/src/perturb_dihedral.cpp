#include "utils/include/serialize.h"
#include "math/include/constants.h"  // PI
#include "math/include/random.h"
#include "math/include/utils_math.h"  // round
#include "monte_carlo/include/perturb_dihedral.h"

namespace feasst {

PerturbDihedral::PerturbDihedral(argtype args)
  : PerturbDihedral(&args) {
  FEASST_CHECK_ALL_USED(args);
}
PerturbDihedral::PerturbDihedral(argtype * args)
  : PerturbDistanceAngle(args) {
  class_name_ = "PerturbDihedral";
}

class MapPerturbDihedral {
 public:
  MapPerturbDihedral() {
    auto obj = MakePerturbDihedral();
    obj->deserialize_map()["PerturbDihedral"] = obj;
  }
};

static MapPerturbDihedral mapper_ = MapPerturbDihedral();

std::shared_ptr<Perturb> PerturbDihedral::create(std::istream& istr) const {
  return std::make_shared<PerturbDihedral>(istr);
}

void PerturbDihedral::precompute(TrialSelect * select, System * system) {
  PerturbDistanceAngle::precompute(select, system);
  ASSERT(select->has_property("dihedral_type"), "cannot obtain dihedral properties");
  dihedral_type_ = feasst::round(select->property("dihedral_type"));
  DEBUG("dihedral_type_ " << dihedral_type_);
}

double PerturbDihedral::random_dihedral_radians(const System& system,
    const TrialSelect * select, Random * random, double * bond_energy) {
  const Dihedral& dihedral = select->configuration(system).unique_type(
    select->particle_type()).dihedral(dihedral_type_);
  ASSERT(dihedral_.deserialize_map().count(dihedral.model()) == 1,
    dihedral.model() << " not found");
  const BondFourBody * model = dihedral_.deserialize_map()[dihedral.model()].get();
  const double beta = system.thermo_params().beta();
  const double radians = model->random_dihedral_radians(dihedral, beta, system.dimension(), random);
  *bond_energy += model->energy(radians, dihedral);
  DEBUG("bond_energy " << *bond_energy);
  return radians;
}

double PerturbDihedral::old_dihedral_energy(const System& system,
    const TrialSelect * select) {
  const Particle& unique_type = select->configuration(system).unique_type(
    select->particle_type());
  ASSERT(dihedral_type_ < unique_type.num_dihedrals(),
    "dihedral_type:" << dihedral_type_ << " >= number of dihedrals:"
    << unique_type.num_dihedrals());
  const Dihedral& dihedral = unique_type.dihedral(dihedral_type_);
  ASSERT(dihedral_.deserialize_map().count(dihedral.model()) == 1,
    dihedral.model() << " not found");
  const BondFourBody * model = dihedral_.deserialize_map()[dihedral.model()].get();
  return model->energy(
    select->mobile().site_positions()[0][0],
    select->anchor_position(0, 0, system),
    select->anchor_position(0, 1, system),
    select->anchor_position(0, 2, system),
    dihedral);
}

void PerturbDihedral::move_once(const bool is_position_held,
    System * system,
    TrialSelect * select,
    Random * random,
    double * bond_energy) {
  DEBUG(class_name());
  if (is_position_held) {
    *bond_energy += old_bond_energy(*system, select);
    *bond_energy += old_angle_energy(*system, select);
    *bond_energy += old_dihedral_energy(*system, select);
  } else {
    const double distance = random_distance(*system, select, random, bond_energy);
    const double angle = random_angle_radians(*system, select, random, bond_energy);
    const double dihedral = random_dihedral_radians(*system, select, random, bond_energy);
    place_dihedral(distance, angle, dihedral, system, select);
  }
  DEBUG("bond_energy dist+ang+dih " << *bond_energy);
}

void PerturbDihedral::place_dihedral(const double distance,
  const double angle,
  const double dihedral,
  System * system,
  TrialSelect * select) {
  DEBUG("placing dihedral at radian angle " << angle);
  const int dimen = select->configuration(*system).dimension();
  ASSERT(dimen == 3, "not implemented for dimen: " << dimen);
  if (origin_.dimension() == 0) origin_.set_to_origin(dimen);
  Select * mobile = select->get_mobile();
  Position * site = mobile->get_site_position(0, 0);
  DEBUG("mobile " << mobile->str());
  DEBUG("old pos " << site->str());

  /*
    For given sites j, k and l (anchors), place site, i, according to bond
    length, angle and dihedral.

    l
     \
      k - j
           \
            i

    r_ij = r_i - r_j points from j to i
    - initialize r_ij as r_jk, normalized and multiplied by bond length
    - define the normal n_2 = r_kj cross r_lk = r_jk cross r_kl
    - rotate r_ij by (180 - angle) about n_2
    - rotate r_ij by (180 - dihedral) about r_jk
   */
  const Position& rj = select->anchor_position(0, 0, *system);
  const Position& rk = select->anchor_position(0, 1, *system);
  const Position& rl = select->anchor_position(0, 2, *system);
  DEBUG("rj: " << rj.str() << " rk: " << rk.str() << " rl: " << rl.str());
  rjk_ = rj;
  rjk_.subtract(rk);
  *site = rjk_;

  // normalize site position, centered at rj, and multiply by bond length.
  site->normalize();
  site->multiply(distance);
  rkl_ = rk;
  rkl_.subtract(rl);
  const Position n1 = rkl_.cross_product(rjk_);
  DEBUG("n1 " << n1.str());
  rot_mat_.axis_angle(n1, radians_to_degrees(PI - angle));
  rot_mat_.rotate(origin_, site);
  rot_mat_.axis_angle(rjk_, radians_to_degrees(dihedral));
  rot_mat_.rotate(origin_, site);
  site->add(rj);  // return to frame of reference
  DEBUG("new pos " << site->str());
  select->get_configuration(system)->update_positions(select->mobile());
}

PerturbDihedral::PerturbDihedral(std::istream& istr)
  : PerturbDistanceAngle(istr) {
  ASSERT(class_name_ == "PerturbDihedral", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(7579 == version, "mismatch version: " << version);
  feasst_deserialize(&dihedral_type_, istr);
}

void PerturbDihedral::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_perturb_dihedral_(ostr);
}

void PerturbDihedral::serialize_perturb_dihedral_(
    std::ostream& ostr) const {
  serialize_perturb_distance_angle_(ostr);
  feasst_serialize_version(7579, ostr);
  feasst_serialize(dihedral_type_, ostr);
}

}  // namespace feasst
