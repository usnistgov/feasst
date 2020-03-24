#include "chain/include/trial_select_segment.h"
#include "utils/include/serialize.h"
#include "math/include/utils_math.h"
#include "math/include/random.h"

namespace feasst {

class MapTrialSelectSegment {
 public:
  MapTrialSelectSegment() {
    auto obj = MakeTrialSelectSegment();
    obj->deserialize_map()["TrialSelectSegment"] = obj;
  }
};

static MapTrialSelectSegment mapper_ = MapTrialSelectSegment();

std::shared_ptr<TrialSelect> TrialSelectSegment::create(std::istream& istr) const {
  return std::make_shared<TrialSelectSegment>(istr);
}

TrialSelectSegment::TrialSelectSegment(std::istream& istr)
  : TrialSelectParticle(istr) {
  // ASSERT(class_name_ == "TrialSelectSegment", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(658 == version, "mismatch version: " << version);
  feasst_deserialize(&max_length_, istr);
}

void TrialSelectSegment::serialize_trial_select_segment_(std::ostream& ostr) const {
  serialize_trial_select_particle_(ostr);
  feasst_serialize_version(658, ostr);
  feasst_serialize(max_length_, ostr);
}

void TrialSelectSegment::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_select_segment_(ostr);
}

TrialSelectSegment::TrialSelectSegment(const argtype& args)
  : TrialSelectParticle(args) {
  class_name_ = "TrialSelectSegment";
  args_.dont_check();
  max_length_ = args_.key("max_length").dflt("-1").integer();
}

bool TrialSelectSegment::random_segment_in_particle(
    const Configuration& config,
    Select * select,
    Random * random,
    const int max_length) {
  random_particle(config, select, random);
  const int num_sites = select->num_sites();
  DEBUG("num_sites " << num_sites);
  if (num_sites <= 1) {
    // HWH note this check prevents error/infinite loop below
    return false;
  }

  // find two unequal sites
  int min = 0;
  int max = min;
  int attempt = 0;
  while (min == max) {
    min = random->uniform(0, num_sites - 1);
    if (max_length == -1) {
      max = random->uniform(0, num_sites - 1);
    } else {
      max = min + random->uniform(-max_length, max_length);
      if (max < 0) {
        max = 0;
      }
      if (max >= num_sites) {
        max = num_sites - 1;
      }
    }
    ++attempt;
    ASSERT(attempt < 1e3, "infinite loop");
  }

  // swap for meaningful min/max
  feasst::sort(&min, &max);

  // remove sites not in min/max, from highest to lowest
  select->remove_last_sites(num_sites - max - 1);
  select->remove_first_sites(min);
  return true;
}

bool TrialSelectSegment::select(const Select& perturbed,
    System* system,
    Random * random) {
  const bool is_found = random_segment_in_particle(
    system->configuration(),
    &mobile_,
    random,
    max_length()
  );
  if (!is_found) {
    return false;
  }
  mobile_original_ = mobile_;
  return true;
}

}  // namespace feasst
