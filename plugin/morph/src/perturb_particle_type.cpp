#include "utils/include/serialize.h"
#include "morph/include/perturb_particle_type.h"

namespace feasst {

PerturbParticleType::PerturbParticleType(argtype args)
  : PerturbParticleType(&args) {
  check_all_used(args);
}
PerturbParticleType::PerturbParticleType(argtype * args) : Perturb(args) {
  class_name_ = "PerturbParticleType";
  disable_tunable_();
  new_particle_type_ = integer("type", args);
}

class MapPerturbParticleType {
 public:
  MapPerturbParticleType() {
    auto obj = MakePerturbParticleType({{"type", "0"}});
    obj->deserialize_map()["PerturbParticleType"] = obj;
  }
};

static MapPerturbParticleType mapper_ = MapPerturbParticleType();

std::shared_ptr<Perturb> PerturbParticleType::create(std::istream& istr) const {
  return std::make_shared<PerturbParticleType>(istr);
}

void PerturbParticleType::perturb(
    System * system,
    TrialSelect * select,
    Random * random,
    const bool is_position_held) {
  if (is_position_held) {
    select->set_trial_state(0);
    return;
  }
  set_particle_type(system, select->mobile(), new_particle_type_);
  set_revert_possible(true, select);
  set_finalize_possible(true, select);
  select->set_trial_state(1);
}

void PerturbParticleType::set_particle_type(
    System * system,
    const Select& select,
    const int type) {
  ASSERT(select.num_particles() == 1, "assumes 1 particle: " << select.str() <<
    " or else haven't implemented revert correctly");
  const int particle_index = select.particle_index(0);
  const Configuration& config = system->configuration();
  old_particle_type_ = config.select_particle(particle_index).type();
  system->get_configuration()->set_particle_type(type, select);
}

void PerturbParticleType::revert(System * system) {
  DEBUG("revert_possible " << revert_possible());
  if (revert_possible()) {
    set_particle_type(system, revert_select()->mobile(), old_particle_type_);
    system->revert(revert_select()->mobile());
  }
}

void PerturbParticleType::finalize(System * system) {
  DEBUG("finalize_possible " << finalize_possible());
  if (finalize_possible()) {
    //system->finalize(finalize_select()->mobile());
  }
}

std::string PerturbParticleType::status_header() const {
  std::stringstream ss;
  return ss.str();
}

std::string PerturbParticleType::status() const {
  std::stringstream ss;
  return ss.str();
}

PerturbParticleType::PerturbParticleType(std::istream& istr)
  : Perturb(istr) {
  ASSERT(class_name_ == "PerturbParticleType", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(4672 == version, "mismatch version: " << version);
  feasst_deserialize(&new_particle_type_, istr);
}

void PerturbParticleType::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_perturb_(ostr);
  feasst_serialize_version(4672, ostr);
  feasst_serialize(new_particle_type_, ostr);
}

}  // namespace feasst
