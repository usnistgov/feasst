#include "chain/include/trial_select_end_segment.h"
#include "utils/include/serialize.h"
#include "math/include/random.h"

namespace feasst {

class MapTrialSelectEndSegment {
 public:
  MapTrialSelectEndSegment() {
    auto obj = MakeTrialSelectEndSegment();
    obj->deserialize_map()["TrialSelectEndSegment"] = obj;
  }
};

static MapTrialSelectEndSegment mapper_ = MapTrialSelectEndSegment();

std::shared_ptr<TrialSelect> TrialSelectEndSegment::create(std::istream& istr) const {
  return std::make_shared<TrialSelectEndSegment>(istr);
}

TrialSelectEndSegment::TrialSelectEndSegment(std::istream& istr)
  : TrialSelectSegment(istr) {
  // ASSERT(class_name_ == "TrialSelectEndSegment", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(456 == version, "mismatch version: " << version);
}

void TrialSelectEndSegment::serialize_trial_select_end_segment_(std::ostream& ostr) const {
  serialize_trial_select_segment_(ostr);
  feasst_serialize_version(456, ostr);
}

void TrialSelectEndSegment::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_select_end_segment_(ostr);
}

void TrialSelectEndSegment::precompute(System * system) {
  TrialSelectSegment::precompute(system);
  anchor_.clear();
  anchor_.add_site(0, 0);
}

bool TrialSelectEndSegment::random_end_segment_in_particle(
    const Configuration& config,
    const int max_length,
    Select * select,
    Random * random,
    bool * is_endpoint_beginning) {
  random_particle(config, select, random);
  const int num_sites = select->num_sites();
  if (num_sites <= 1) {
    DEBUG("num sites(" << num_sites << ") not large enough");
    return false;
  }

  // select a random site
  int site = -1;
  if (max_length == -1) {
    site = random->uniform(0, num_sites - 1);

    DEBUG("site " << site << " num " << num_sites);

    // randomly decide which endpoint to keep in selection
    if (site == 0) {
      *is_endpoint_beginning = false;
    } else if (site == num_sites - 1) {
      *is_endpoint_beginning = true;
    } else {
      if (random->coin_flip()) {
        *is_endpoint_beginning = false;
      } else {
        *is_endpoint_beginning = true;
      }
    }
  } else {
    ASSERT(max_length > 0, "max_length(" << max_length <<") should be >0 "
      << "or no segment will be selected");
    if (random->coin_flip()) {
      *is_endpoint_beginning = false;
      site = random->uniform(num_sites - max_length, num_sites - 1);
    } else {
      *is_endpoint_beginning = true;
      site = random->uniform(0, max_length - 1);
    }
  }

  DEBUG("beginning? " << *is_endpoint_beginning);
  if (*is_endpoint_beginning) {
    select->remove_last_sites(num_sites - site - 1);
  } else {
    select->remove_first_sites(site);
  }
  DEBUG("num " << num_sites << " indices " << select->str());
  return true;
}

bool TrialSelectEndSegment::select(const Select& perturbed,
    System* system,
    Random * random) {
  bool is_endpoint_beginning;
  const bool is_found = random_end_segment_in_particle(
    system->configuration(),
    max_length(),
    &mobile_,
    random,
    &is_endpoint_beginning
  );
  if (!is_found) {
    return false;
  }
  update_anchor(is_endpoint_beginning, system);
  mobile_original_ = mobile_;
  return true;
}

void TrialSelectEndSegment::update_anchor(const bool is_endpoint_beginning,
  const System * system) {
  int select_index = -1;
  if (is_endpoint_beginning) {
    select_index = mobile_.num_sites() - 1;
  } else {
    select_index = 0;
  }
  DEBUG("is_endpoint_beginning, " << is_endpoint_beginning);
  DEBUG("site index " << select_index);
  anchor_.set_site(0, 0, mobile_.site_indices()[0][select_index]);
  anchor_.set_particle(0, mobile_.particle_indices()[0]);
}

}  // namespace feasst
