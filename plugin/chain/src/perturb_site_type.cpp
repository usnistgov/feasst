#include "utils/include/serialize.h"
#include "utils/include/arguments.h"
#include "configuration/include/select.h"
#include "configuration/include/particle_factory.h"
#include "configuration/include/configuration.h"
#include "system/include/system.h"
#include "monte_carlo/include/trial_select.h"
#include "chain/include/perturb_site_type.h"

namespace feasst {

PerturbSiteType::PerturbSiteType(argtype args) : PerturbSiteType(&args) {
  feasst_check_all_used(args);
}
PerturbSiteType::PerturbSiteType(argtype * args) : Perturb(args) {
  class_name_ = "PerturbSiteType";
  disable_tunable_();
  new_site_type_ = integer("type", args);
}

class MapPerturbSiteType {
 public:
  MapPerturbSiteType() {
    auto obj = MakePerturbSiteType({{"type", "0"}});
    obj->deserialize_map()["PerturbSiteType"] = obj;
  }
};

static MapPerturbSiteType mapper_ = MapPerturbSiteType();

std::shared_ptr<Perturb> PerturbSiteType::create(std::istream& istr) const {
  return std::make_shared<PerturbSiteType>(istr);
}

void PerturbSiteType::perturb(
    System * system,
    TrialSelect * select,
    Random * random,
    const bool is_position_held,
    Acceptance * acceptance) {
  if (is_position_held) {
    select->set_trial_state(0);
    return;
  }
  set_site_type(system, *select, new_site_type_);
  set_revert_possible(true, select);
  set_finalize_possible(true, select);
  select->set_trial_state(1);
}

void PerturbSiteType::set_site_type(
    System * system,
    const TrialSelect& tselect,
    const int type) {
  const Select& select = tselect.mobile();
  ASSERT(select.num_particles() == 1, "assumes 1 particle: "
    << select.num_particles());
  ASSERT(select.num_sites() == 1, "assumes 1 site: " << select.num_sites());
  const Configuration& config = tselect.configuration(*system);
  const Particle& part = config.select_particle(select.particle_index(0));
  old_site_type_ = part.site(select.site_index(0, 0)).type();
  old_particle_type_ = part.type();
  const int site = select.site_index(0, 0);
  tselect.get_configuration(system)->set_site_type(old_particle_type_, site, type);
}

void PerturbSiteType::revert(System * system) {
  DEBUG("revert_possible " << revert_possible());
  if (revert_possible()) {
    set_site_type(system, *revert_select(), old_site_type_);
    system->revert(revert_select()->mobile());
  }
}

void PerturbSiteType::finalize(System * system) {
  DEBUG("finalize_possible " << finalize_possible());
  if (finalize_possible()) {
//    system->finalize(finalize_select()->mobile());
  }
}

std::string PerturbSiteType::status_header() const {
  std::stringstream ss;
  return ss.str();
}

std::string PerturbSiteType::status() const {
  std::stringstream ss;
  return ss.str();
}

PerturbSiteType::PerturbSiteType(std::istream& istr)
  : Perturb(istr) {
  ASSERT(class_name_ == "PerturbSiteType", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(953 == version, "mismatch version: " << version);
  feasst_deserialize(&new_site_type_, istr);
}

void PerturbSiteType::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_perturb_(ostr);
  feasst_serialize_version(953, ostr);
  feasst_serialize(new_site_type_, ostr);
}

void PerturbSiteType::precompute(TrialSelect * select, System * system) {
  select->set_ghost(true);
}

}  // namespace feasst
