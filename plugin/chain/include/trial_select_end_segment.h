
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
  void precompute(System * system) override;

  /// Select all sites between a random endpoint and a randomly selectioned site in a randomly selected particle in group.
  /// Return true if the endpoint is at the beginning.
  bool random_end_segment_in_particle(const Configuration& config,
      const int max_length,
      Select * select,
      Random * random,
      /// Set the maximum length of the segment.
      /// If -1 (default), consider all possible lengths.
      bool * is_endpoint_beginning);

  bool select(const Select& perturbed, System* system, Random * random) override;

  virtual void update_anchor(const bool is_endpoint_beginning,
    const System * system);

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
