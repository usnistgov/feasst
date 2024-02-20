#include "utils/include/serialize.h"
#include "morph/include/perturb_particle_type.h"

namespace feasst {

PerturbParticleType::PerturbParticleType(argtype args)
  : PerturbParticleType(&args) {
  FEASST_CHECK_ALL_USED(args);
}
PerturbParticleType::PerturbParticleType(argtype * args) : Perturb(args) {
  class_name_ = "PerturbParticleType";
  rotate_.set_tunable_min_and_max(-1 - NEAR_ZERO, -1 + NEAR_ZERO);
  rotate_.set_tunable(-1);
  rotate_.disable_tunable_();
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
    const bool is_position_held,
    Acceptance * acceptance) {
  if (is_position_held) {
    select->set_trial_state(0);
    return;
  }
  set_particle_type(system, select->mobile(), new_particle_type_);
  const int num_sites = select->mobile().num_sites();
  if (num_sites > 1) {
    // use particle_type to set the coordinates, and randomly rotate them about the anchor.
    // then set to position of existing anchor
    const Configuration& config = system->configuration();
    const Particle& ref_part = config.particle_type(new_particle_type_);
    ASSERT(num_sites == ref_part.num_sites(), "num sites mismatch");
    ASSERT(select->mobile().num_particles() == 1,
      "Implemented for only one particle if multiple sites.")
    const Position& ref_site_pos_old = select->mobile().site_positions()[0][0];
    //const int part_index = select->mobile().particle_indices()[0];
    const Position& ref_site_pos_new = ref_part.site(0).position();
    //if (tmp_pos_.size() == 0) {
    //  tmp_pos_.set_as_origin(ref_site_pos_old.dimension());
    //}
    for (int site = 1; site < num_sites; ++site) {
      tmp_pos_ = ref_part.site(site).position();
      tmp_pos_.subtract(ref_site_pos_new);
      tmp_pos_.add(ref_site_pos_old);
      select->get_mobile()->set_site_position(0, site, tmp_pos_);
    }
    rotate_.move(is_position_held, system, select, random, acceptance);
  }
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
    Configuration* config = system->get_configuration();
    config->update_positions(revert_select()->mobile_original(),
      // don't wrap if reverting
      false);
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
  feasst_deserialize_fstobj(&rotate_, istr);
}

void PerturbParticleType::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_perturb_(ostr);
  feasst_serialize_version(4672, ostr);
  feasst_serialize(new_particle_type_, ostr);
  feasst_serialize_fstobj(rotate_, ostr);
}

}  // namespace feasst
