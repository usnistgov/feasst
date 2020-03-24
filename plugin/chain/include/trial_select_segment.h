
#ifndef FEASST_CHAIN_TRIAL_SELECT_SEGMENT_H_
#define FEASST_CHAIN_TRIAL_SELECT_SEGMENT_H_

#include "monte_carlo/include/trial_select_particle.h"

namespace feasst {

/// Select a random segment.
// HWH optimzie by settting endpoints as anchors.
class TrialSelectSegment : public TrialSelectParticle {
 public:
  /**
    args:
    - max_length : maximum length of selected segment. If -1 (default), then
      randomly select all possible lengths.
   */
  TrialSelectSegment(const argtype& args = argtype());

  /// Return the maximum length.
  int max_length() const { return max_length_; }

  /// Select all sites between two randomly selected sites in a randomly selected particle in group.
  /// Return true if selection is valid
  bool random_segment_in_particle(
      const Configuration& config,
      Select * select,
      Random * random,
      /// Set the maximum length of the segment.
      /// If -1 (default), consider all possible lengths.
      const int max_length = -1);

  bool select(const Select& perturbed,
    System* system,
    Random * random) override;

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
