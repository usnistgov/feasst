
#ifndef FEASST_CHAIN_SELECT_END_SEGMENT_H_
#define FEASST_CHAIN_SELECT_END_SEGMENT_H_

#include "chain/include/select_segment.h"

namespace feasst {

// HWH optimize, set anchor one site next from selection.
/// Select an end segment.
/// Set the anchor as the other end of the selection from the end point.
class SelectEndSegment : public SelectSegment {
 public:
  //@{
  /** @name Arguments
   */

  /**
    args:
    - SelectSegment arguments.
   */
  explicit SelectEndSegment(argtype args = argtype());
  explicit SelectEndSegment(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

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
  explicit SelectEndSegment(std::istream& istr);
  virtual ~SelectEndSegment() {}

  //@}
 protected:
  void serialize_select_end_segment_(std::ostream& ostr) const;
};

inline std::shared_ptr<SelectEndSegment> MakeSelectEndSegment(
    const argtype &args = argtype()) {
  return std::make_shared<SelectEndSegment>(args);
}

}  // namespace feasst

#endif  // FEASST_CHAIN_SELECT_END_SEGMENT_H_
