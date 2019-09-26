
#ifndef FEASST_CHAIN_TRIAL_SELECT_SEGMENT_H_
#define FEASST_CHAIN_TRIAL_SELECT_SEGMENT_H_

#include "monte_carlo/include/trial_select_particle.h"

namespace feasst {

/// Select a random segment.
// HWH optimzie by settting endpoints as anchors.
class TrialSelectSegment : public TrialSelectParticle {
 public:
  TrialSelectSegment(
    /**
      max_length : maximum length of selected segment. If -1 (default), then
        randomly select all possible lengths.
     */
    const argtype& args = argtype()) : TrialSelectParticle(args) {
    class_name_ = "TrialSelectSegment";
    Arguments args_(args);
    args_.dont_check();
    max_length_ = args_.key("max_length").dflt("-1").integer();
  }

  int max_length() const { return max_length_; }

  /// Select all sites between two randomly selected sites in a randomly selected particle in group.
  void random_segment_in_particle(
      const Configuration& config,
      SelectPosition * select,
      Random * random,
      /// Set the maximum length of the segment.
      /// If -1 (default), consider all possible lengths.
      const int max_length = -1
    ) {
    random_particle(config, select, random);
    const int num_sites = select->num_sites();
    if (num_sites <= 1) {
      return; // HWH note this check prevents error/infinite loop below
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
    sort(&min, &max);

    // remove sites not in min/max, from highest to lowest
    select->remove_last_sites(num_sites - max - 1);
    select->remove_first_sites(min);
  }

  bool select(const Select& perturbed, System* system, Random * random) override {
    random_segment_in_particle(
      system->configuration(),
      &mobile_,
      random,
      max_length()
    );
    mobile_original_ = mobile_;
    return true;
  }

  std::shared_ptr<TrialSelect> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit TrialSelectSegment(std::istream& istr);
  virtual ~TrialSelectSegment() {}

 protected:
  void serialize_trial_select_segment_(std::ostream& ostr) const;

 private:
  int max_length_;
};

inline std::shared_ptr<TrialSelectSegment> MakeTrialSelectSegment(
    const argtype &args = argtype()) {
  return std::make_shared<TrialSelectSegment>(args);
}

}  // namespace feasst

#endif  // FEASST_CHAIN_TRIAL_SELECT_SEGMENT_H_
