#include "chain/include/trial_select_site_of_type.h"
#include "utils/include/serialize.h"
#include "math/include/random.h"

namespace feasst {

class MapTrialSelectSiteOfType {
 public:
  MapTrialSelectSiteOfType() {
    auto obj = MakeTrialSelectSiteOfType({{"site_type", "0"}});
    obj->deserialize_map()["TrialSelectSiteOfType"] = obj;
  }
};

static MapTrialSelectSiteOfType mapper_ = MapTrialSelectSiteOfType();

std::shared_ptr<TrialSelect> TrialSelectSiteOfType::create(std::istream& istr) const {
  return std::make_shared<TrialSelectSiteOfType>(istr);
}

TrialSelectSiteOfType::TrialSelectSiteOfType(std::istream& istr)
  : TrialSelect(istr) {
  // ASSERT(class_name_ == "TrialSelectSiteOfType", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(272 == version, "mismatch version: " << version);
  feasst_deserialize(&site_type_, istr);
}

void TrialSelectSiteOfType::serialize_trial_select_segment_(std::ostream& ostr) const {
  serialize_trial_select_(ostr);
  feasst_serialize_version(272, ostr);
  feasst_serialize(site_type_, ostr);
}

void TrialSelectSiteOfType::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_select_segment_(ostr);
}

TrialSelectSiteOfType::TrialSelectSiteOfType(const argtype& args)
  : TrialSelect(args) {
  class_name_ = "TrialSelectSiteOfType";
  args_.dont_check();
  site_type_ = args_.key("site_type").integer();
}

int TrialSelectSiteOfType::random_site_in_particle(
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

bool TrialSelectSiteOfType::select(const Select& perturbed,
    System* system,
    Random * random) {
  const int num = random_site_in_particle(
    system->configuration(),
    &mobile_,
    random);
  if (num <= 0) return false;
  set_probability(1./static_cast<double>(num));
  mobile_original_ = mobile_;
  return true;
}

}  // namespace feasst
