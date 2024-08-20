#include "utils/include/serialize.h"
#include "utils/include/arguments.h"
#include "configuration/include/select.h"
#include "configuration/include/particle_factory.h"
#include "configuration/include/configuration.h"
#include "monte_carlo/include/trial_select.h"
#include "chain/include/perturb_distance_angle_connector.h"

namespace feasst {

PerturbDistanceAngleConnector::PerturbDistanceAngleConnector(argtype args)
  : PerturbDistanceAngleConnector(&args) {
  feasst_check_all_used(args);
}
PerturbDistanceAngleConnector::PerturbDistanceAngleConnector(argtype * args)
  : PerturbDistanceAngle(args) {
  class_name_ = "PerturbDistanceAngleConnector";
}

class MapPerturbDistanceAngleConnector {
 public:
  MapPerturbDistanceAngleConnector() {
    auto obj = MakePerturbDistanceAngleConnector();
    obj->deserialize_map()["PerturbDistanceAngleConnector"] = obj;
  }
};

static MapPerturbDistanceAngleConnector mapper_ = MapPerturbDistanceAngleConnector();

std::shared_ptr<Perturb> PerturbDistanceAngleConnector::create(std::istream& istr) const {
  return std::make_shared<PerturbDistanceAngleConnector>(istr);
}

void PerturbDistanceAngleConnector::move_once(const bool is_position_held,
    System * system,
    TrialSelect * select,
    Random * random,
    double * bond_energy) {
  PerturbDistanceAngle::move_once(is_position_held, system, select, random, bond_energy);
  if (is_position_held) return;
  // find vector between mobile (anisotropic center) and anchor (connector)
  // find this vector for both internal (fstprt) and external frame of reference (mobile)
  const Configuration& config = select->configuration(*system);
  const Particle& part = config.select_particle(select->mobile().particle_index(0));
  const Particle& part_type = config.particle_type(part.type());
  Position internal = part_type.site(select->anchor().site_index(0, 0)).position();
  internal.subtract(part_type.site(select->mobile().site_index(0, 0)).position());
  DEBUG("internal " << internal.str());
  Position external = select->anchor_position(0, 0, *system);
  external.subtract(select->mobile().site_positions()[0][0]);
  DEBUG("external " << external.str());

  // obtain rotation matrix which converts internal to current
  rot_mat_.vector_onto_vector(internal, external, &tmp_pos_, &tmp_mat_);

  // obtain euler angles from rotation matrix
  euler_.set(rot_mat_);

  // update the euler angles in mobile selection
  select->get_mobile()->set_euler(0, 0, euler_);

  // update the euler angles in configuration?
  select->get_configuration(system)->update_positions(select->mobile());
  // optimization is possible here, as update_positions is called twice.
  // once in base class, and now again here.
}

PerturbDistanceAngleConnector::PerturbDistanceAngleConnector(std::istream& istr)
  : PerturbDistanceAngle(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(7548 == version, "mismatch version: " << version);
}

void PerturbDistanceAngleConnector::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_perturb_distance_angle_connector_(ostr);
}

void PerturbDistanceAngleConnector::serialize_perturb_distance_angle_connector_(
    std::ostream& ostr) const {
  serialize_perturb_distance_angle_(ostr);
  feasst_serialize_version(7548, ostr);
}

}  // namespace feasst
