#include <cmath>
#include "utils/include/serialize.h"
#include "math/include/constants.h"
#include "math/include/utils_math.h"
#include "math/include/random.h"
#include "monte_carlo/include/perturb_distance.h"

namespace feasst {

PerturbDistance::PerturbDistance(argtype args) : PerturbDistance(&args) {
  FEASST_CHECK_ALL_USED(args);
}
PerturbDistance::PerturbDistance(argtype * args) : PerturbMove(args) {
  class_name_ = "PerturbDistance";
  if (!boolean("enable_tunable", args, false)) {
    disable_tunable_();
  } else {
    set_tunable_min_and_max(2*NEAR_ZERO, 180.*(1 + NEAR_ZERO));
  }
  potential_acceptance_ = integer("potential_acceptance", args, -1);
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
  ASSERT(select->has_property("bond_type"), "cannot obtain bond properties");
  bond_type_ = feasst::round(select->property("bond_type"));
}

void PerturbDistance::move(const bool is_position_held,
                           System * system,
                           TrialSelect * select,
                           Random * random) {
  double bond_energy = 0.;
  DEBUG("potential_acceptance_ " << potential_acceptance_);
  if (potential_acceptance_ == -1) {
    move_once(is_position_held, system, select, random, &bond_energy);
    select->add_exclude_energy(bond_energy);
    return;
  }
  int max_attempt = 1e6;
  const double beta = system->thermo_params().beta();
  for (int attempt = 0; attempt < max_attempt; ++attempt) {
    move_once(is_position_held, system, select, random, &bond_energy);
    Potential * poten = system->get_potential(potential_acceptance_);
    DEBUG("sel " << select->mobile().str());
    const double energy = poten->select_energy(select->mobile(),
      select->get_configuration(system));
    DEBUG("energy " << energy);
    if (is_position_held || random->uniform() < std::exp(-beta*energy)) {
      select->add_exclude_energy(bond_energy);
      DEBUG("accepted");
      return;
    }
  }
  FATAL("max_attempt: " << max_attempt << " reached");
}

double PerturbDistance::old_bond_energy(const System& system,
    const TrialSelect * select) {
  // the following five lines were copy-pasted from random_distance
  const Bond& bond = select->configuration(system).unique_types().particle(
    select->particle_type()).bond(bond_type_);
  ASSERT(bond_.deserialize_map().count(bond.model()) == 1,
    bond.model() << " not found");
  const BondTwoBody * model = bond_.deserialize_map()[bond.model()].get();
  return model->energy(
    select->mobile().site_positions()[0][0],
    select->anchor_position(0, 0, system),
    bond);
}

void PerturbDistance::move_once(const bool is_position_held,
    System * system,
    TrialSelect * select,
    Random * random,
    double * bond_energy) {
  DEBUG(class_name());
  if (is_position_held) {
    *bond_energy = old_bond_energy(*system, select);
    return;
  }
  Position * site = select->get_mobile()->get_site_position(0, 0);
  const Position& anchor_pos = select->anchor_position(0, 0, *system);
  DEBUG("mobile " << select->mobile().str());
  DEBUG("old pos " << site->str());
  if (tunable().is_enabled()) {
    //INFO("tunable enabled? " << tunable().is_enabled());
    //INFO("tunable value? " << tunable().value());
    const double tune = tunable().value();
    ASSERT(std::abs(tune) > NEAR_ZERO, "max angle is too small");
    random->rotation(site->dimension(), &axis_tmp_, &rot_mat_tmp_, tune);
    site->subtract(anchor_pos);
    site->normalize();
    site->multiply(random_distance(*system, select, random, bond_energy));
    DEBUG("final dist pert bond_energy " << *bond_energy);
    if (origin_tmp_.dimension() == 0) {
      origin_tmp_.set_to_origin(site->dimension());
    }
    rot_mat_tmp_.rotate(origin_tmp_, site);
  } else {
    random->unit_sphere_surface(site);
    site->multiply(random_distance(*system, select, random, bond_energy));
    DEBUG("final dist pert bond_energy " << *bond_energy);
  }
  site->add(anchor_pos);
  DEBUG("new pos " << site->str());
  select->get_configuration(system)->update_positions(select->mobile());
}

PerturbDistance::PerturbDistance(std::istream& istr)
  : PerturbMove(istr) {
  // HWH can't check this if this is a base class
  // ASSERT(class_name_ == "PerturbDistance", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(228 == version, "mismatch version: " << version);
  feasst_deserialize(&bond_type_, istr);
  feasst_deserialize(&potential_acceptance_, istr);
}

void PerturbDistance::serialize_perturb_distance_(std::ostream& ostr) const {
  serialize_perturb_(ostr);
  feasst_serialize_version(228, ostr);
  feasst_serialize(bond_type_, ostr);
  feasst_serialize(potential_acceptance_, ostr);
}

void PerturbDistance::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_perturb_distance_(ostr);
}

double PerturbDistance::random_distance(const System& system,
    const TrialSelect* select,
    Random * random,
    double * bond_energy) {
  const Bond& bond = select->configuration(system).unique_types().particle(
    select->particle_type()).bond(bond_type_);
  ASSERT(bond_.deserialize_map().count(bond.model()) == 1,
    bond.model() << " not found");
  const BondTwoBody * model = bond_.deserialize_map()[bond.model()].get();
  const double beta = system.thermo_params().beta();
  const double dist = model->random_distance(bond, beta, system.dimension(), random);
  *bond_energy += model->energy(dist, bond);
  DEBUG("bond_energy " << *bond_energy);
  return dist;
}

}  // namespace feasst
