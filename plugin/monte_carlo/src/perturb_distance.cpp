#include <cmath>
#include "utils/include/serialize.h"
#include "math/include/constants.h"
#include "math/include/utils_math.h"
#include "math/include/random.h"
#include "monte_carlo/include/perturb_distance.h"

namespace feasst {

PerturbDistance::PerturbDistance(argtype args) : PerturbDistance(&args) {
  check_all_used(args);
}
PerturbDistance::PerturbDistance(argtype * args) : PerturbMove(args) {
  class_name_ = "PerturbDistance";
  disable_tunable_();
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
  const int bond_type = feasst::round(select->property("bond_type"));
  const Bond& bond = system->configuration().unique_types().particle(
    select->particle_type()).bond(bond_type);
  if (bond.has_property("length")) {
    distance_ = bond.property("length");
  }
  if (bond.has_property("spring_constant")) {
    spring_constant_ = bond.property("spring_constant");
  }
  if (bond.has_property("maximum_length")) {
    maximum_length_ = bond.property("maximum_length");
  }
  if (bond.has_property("exponent")) {
    exponent_ = bond.property("exponent");
  }
}

void PerturbDistance::move(System * system,
                           TrialSelect * select,
                           Random * random) {
  if (potential_acceptance_ == -1) {
    move_once_(system, select, random);
    return;
  }
  int max_attempt = 1e6;
  const double beta = system->thermo_params().beta();
  for (int attempt = 0; attempt < max_attempt; ++attempt) {
    move_once_(system, select, random);
    const double energy = system->get_potential(potential_acceptance_)->select_energy(select->mobile(), system->get_configuration());
    if (random->uniform() < std::exp(-beta*energy)) {
      return;
    }
  }
  FATAL("max_attempt: " << max_attempt << " reached");
}

void PerturbDistance::move_once_(System * system,
    TrialSelect * select,
    Random * random) {
  DEBUG(class_name());
  Select * mobile = select->get_mobile();
  Position * site = mobile->get_site_position(0, 0);
  DEBUG("mobile " << mobile->str());
  DEBUG("old pos " << site->str());
  random->unit_sphere_surface(site);
  site->multiply(random_distance(random,
    system->thermo_params().beta(),
    system->dimension()));
  site->add(select->anchor_position(0, 0, *system));
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
  feasst_deserialize(&spring_constant_, istr);
  feasst_deserialize(&maximum_length_, istr);
  feasst_deserialize(&exponent_, istr);
  feasst_deserialize(&potential_acceptance_, istr);
}

void PerturbDistance::serialize_perturb_distance_(std::ostream& ostr) const {
  serialize_perturb_(ostr);
  feasst_serialize_version(228, ostr);
  feasst_serialize(distance_, ostr);
  feasst_serialize(spring_constant_, ostr);
  feasst_serialize(maximum_length_, ostr);
  feasst_serialize(exponent_, ostr);
  feasst_serialize(potential_acceptance_, ostr);
}

void PerturbDistance::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_perturb_distance_(ostr);
}

double PerturbDistance::random_distance(Random * random,
    const double beta,
    const int dimension) const {
  if (std::abs(spring_constant_ + 1) < NEAR_ZERO) {
    return distance_;
  }
  double spring = spring_constant_/beta;
  if (std::abs(spring_constant_) < NEAR_ZERO) spring = 0.;
  double max_length = maximum_length_;
  if (std::abs(max_length + 1) < NEAR_ZERO) {
    max_length = 2*distance_;
  }
  double expo = exponent_;
  if (std::abs(expo + 1) < NEAR_ZERO) {
    expo = 2;
  }
  return random->bond_length(distance_, max_length, spring, expo, dimension);
}

}  // namespace feasst
