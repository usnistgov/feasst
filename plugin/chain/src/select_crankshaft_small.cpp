#include "utils/include/serialize.h"
#include "math/include/utils_math.h"
#include "math/include/random.h"
#include "chain/include/select_crankshaft_small.h"

namespace feasst {

class MapSelectCrankshaftSmall {
 public:
  MapSelectCrankshaftSmall() {
    auto obj = MakeSelectCrankshaftSmall({{"site", "0"},
                                          {"anchor_site0", "1"},
                                          {"anchor_site1", "2"}});
    obj->deserialize_map()["SelectCrankshaftSmall"] = obj;
  }
};

static MapSelectCrankshaftSmall mapper_ = MapSelectCrankshaftSmall();

std::shared_ptr<TrialSelect> SelectCrankshaftSmall::create(std::istream& istr) const {
  return std::make_shared<SelectCrankshaftSmall>(istr);
}

SelectCrankshaftSmall::SelectCrankshaftSmall(std::istream& istr)
  : TrialSelectParticle(istr) {
  // ASSERT(class_name_ == "SelectCrankshaftSmall", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(9287 == version, "mismatch version: " << version);
}

void SelectCrankshaftSmall::serialize_select_crankshaft_small_(std::ostream& ostr) const {
  serialize_trial_select_particle_(ostr);
  feasst_serialize_version(9287, ostr);
}

void SelectCrankshaftSmall::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_select_crankshaft_small_(ostr);
}

SelectCrankshaftSmall::SelectCrankshaftSmall(argtype args) : SelectCrankshaftSmall(&args) {
  FEASST_CHECK_ALL_USED(args);
}
SelectCrankshaftSmall::SelectCrankshaftSmall(argtype * args) : TrialSelectParticle(args) {
  class_name_ = "SelectCrankshaftSmall";
  ASSERT(mobile().num_sites() == 1, "mobile().num_sites() == " <<
    mobile().num_sites() << ". site argument required");
  DEBUG("parse mobile sites");
  std::string start = "site";
  std::stringstream key;
  int num = 1;
  key << start << num;
  while (used(key.str(), *args)) {
    get_mobile()->add_site(0, integer(key.str(), args));
    ++num;
    ASSERT(num < 1e8, "num(" << num << ") is very high. Infinite loop?");
    key.str("");
    key << start << num;
  }
  get_anchor()->clear();
  get_anchor()->add_site(0, integer("anchor_site0", args));
  get_anchor()->add_site(0, integer("anchor_site1", args));
}

bool SelectCrankshaftSmall::select(const Select& perturbed,
    System* system,
    Random * random) {
  Configuration * config = system->get_configuration();
  const int group_index = config->particle_type_to_group_create(particle_type());
  const int num = config->num_particles(group_index);
  if (num <= 0) return false;
  const int index = random->uniform(0, num - 1);
  const Select& select = config->group_select(group_index);
  const int particle_index = select.particle_index(index);
  set_probability_(1./static_cast<double>(num));
  mobile_.set_particle(0, particle_index);
  anchor_.set_particle(0, particle_index);
  mobile_.load_positions(config->particles());
  mobile_original_ = mobile_;
  return true;
}

}  // namespace feasst
