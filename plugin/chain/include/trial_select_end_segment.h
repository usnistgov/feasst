
#ifndef FEASST_CHAIN_TRIAL_SELECT_END_SEGMENT_H_
#define FEASST_CHAIN_TRIAL_SELECT_END_SEGMENT_H_

#include "chain/include/trial_select_segment.h"

namespace feasst {

/// Select an end segment.
/// Set the anchor as the other end of the selection from the end point.
// HWH optimize, set anchor one site next from selection.
class TrialSelectEndSegment : public TrialSelectSegment {
 public:
  TrialSelectEndSegment(const argtype& args = argtype()) : TrialSelectSegment(args) {
    class_name_ = "TrialSelectEndSegment";
  }
  void precompute(System * system) override {
    anchor_.clear();
    anchor_.add_site(0, 0);
  }

  /// Select all sites between a random endpoint and a randomly selectioned site in a randomly selected particle in group.
  /// Return true if the endpoint is at the beginning.
  bool random_end_segment_in_particle(const Configuration& config,
      SelectPosition * select,
      Random * random,
      /// Set the maximum length of the segment.
      /// If -1 (default), consider all possible lengths.
      const int max_length = -1
      ) {
    random_particle(config, select, random);
    const int num_sites = select->num_sites();
    // HWH note this check prevents error/infinite loop below
    if (num_sites <= 1) {
      DEBUG("num sites(" << num_sites << ") not large enough");
      return false;
    }

    // select a random site
    int site = -1;
    bool is_endpoint_beginning;
    if (max_length == -1) {
      site = random->uniform(0, num_sites - 1);

      DEBUG("site " << site << " num " << num_sites);

      // randomly decide which endpoint to keep in selection
      if (site == 0) {
        is_endpoint_beginning = false;
      } else if (site == num_sites - 1) {
        is_endpoint_beginning = true;
      } else {
        if (random->coin_flip()) {
          is_endpoint_beginning = false;
        } else {
          is_endpoint_beginning = true;
        }
      }
    } else {
      ASSERT(max_length > 0, "max_length(" << max_length <<") should be >0 "
        << "or no segment will be selected");
      if (random->coin_flip()) {
        is_endpoint_beginning = false;
        site = random->uniform(num_sites - max_length, num_sites - 1);
      } else {
        is_endpoint_beginning = true;
        site = random->uniform(0, max_length - 1);
      }
    }

    DEBUG("beginning? " << is_endpoint_beginning);
    if (is_endpoint_beginning) {
      select->remove_last_sites(num_sites - site - 1);
    } else {
      select->remove_first_sites(site);
    }
    DEBUG("num " << num_sites << " indices " << select->str());
    return is_endpoint_beginning;
  }

  bool select(const Select& perturbed, System* system, Random * random) override {
    bool is_endpoint_beginning = random_end_segment_in_particle(
      system->configuration(),
      &mobile_,
      random,
      max_length()
    );
    if (mobile_.num_sites() <= 0) {
      return false;
    }
    update_anchor(is_endpoint_beginning, system);
    mobile_original_ = mobile_;
    return true;
  }

  virtual void update_anchor(const bool is_endpoint_beginning,
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

  std::shared_ptr<TrialSelect> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit TrialSelectEndSegment(std::istream& istr);
  virtual ~TrialSelectEndSegment() {}

 protected:
  void serialize_trial_select_end_segment_(std::ostream& ostr) const;
};

inline std::shared_ptr<TrialSelectEndSegment> MakeTrialSelectEndSegment(
    const argtype &args = argtype()) {
  return std::make_shared<TrialSelectEndSegment>(args);
}

}  // namespace feasst

#endif  // FEASST_CHAIN_TRIAL_SELECT_END_SEGMENT_H_
