#include "utils/include/arguments.h"
#include "utils/include/serialize.h"
#include "math/include/random.h"
#include "configuration/include/select.h"
#include "configuration/include/particle_factory.h"
#include "configuration/include/configuration.h"
#include "chain/include/select_site_of_type.h"

namespace feasst {

FEASST_MAPPER(SelectSiteOfType, argtype({{"site_type", "0"}}));

std::shared_ptr<TrialSelect> SelectSiteOfType::create(std::istream& istr) const {
  return std::make_shared<SelectSiteOfType>(istr);
}

SelectSiteOfType::SelectSiteOfType(std::istream& istr)
  : TrialSelect(istr) {
  // ASSERT(class_name_ == "SelectSiteOfType", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(272 == version, "mismatch version: " << version);
  feasst_deserialize(&site_type_, istr);
}

void SelectSiteOfType::serialize_select_segment_(std::ostream& ostr) const {
  serialize_trial_select_(ostr);
  feasst_serialize_version(272, ostr);
  feasst_serialize(site_type_, ostr);
}

void SelectSiteOfType::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_select_segment_(ostr);
}

SelectSiteOfType::SelectSiteOfType(argtype args) : SelectSiteOfType(&args) {
  feasst_check_all_used(args);
}
SelectSiteOfType::SelectSiteOfType(argtype * args) : TrialSelect(args) {
  class_name_ = "SelectSiteOfType";
  site_type_ = integer("site_type", args);
}

int SelectSiteOfType::random_site_in_particle(
    const Configuration& config,
    Select * select,
    Random * random) {
  ASSERT(config.num_site_types() > site_type(),
    "site_type: " << site_type() << " is not present in system.");
  DEBUG("group_index " << group_index());
  const Select& group = config.group_select(group_index());
  const int num_particles = group.num_particles();
  if (num_particles == 0) {
    DEBUG("no particles");
    return 0;
  }
  const int pindex = random->uniform(0, num_particles - 1);
  const int particle_index = group.particle_index(pindex);
  DEBUG("particle_index " << particle_index);
  const Particle& part = config.select_particle(particle_index);
  const int num_sites_of_type = part.num_sites_of_type(site_type());
  if (num_sites_of_type == 0) {
    DEBUG("no sites of type");
    return 0;
  }
  const int sindex = random->uniform(0, num_sites_of_type - 1);

  // loop through each site in particle to select the one of type
  int count = -1;
  int site_index = -1;
  bool terminate = false;
  while (!terminate) {
    ++site_index;
    if (part.site(site_index).type() == site_type()) {
      ++count;
      if (count == sindex) {
        terminate = true;
      }
    }
    ASSERT(site_index < part.num_sites(), "infinite loop");
  }
  if (select->num_particles() != 1 || select->num_sites(0) != 1) {
    select->clear();
    select->add_site(0, 0);
  }
  select->set_particle(0, particle_index);
  select->set_site(0, 0, site_index);
  DEBUG("selected: " << select->str());
  return num_sites_of_type;
}

bool SelectSiteOfType::select(const Select& perturbed,
    System* system,
    Random * random) {
  const int num = random_site_in_particle(
    configuration(*system),
    get_mobile(),
    random);
  if (num <= 0) return false;
  set_probability_(1./static_cast<double>(num));
  set_mobile_original(system);
  return true;
}

}  // namespace feasst
