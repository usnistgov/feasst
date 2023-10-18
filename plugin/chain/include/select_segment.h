
#ifndef FEASST_CHAIN_SELECT_SEGMENT_H_
#define FEASST_CHAIN_SELECT_SEGMENT_H_

#include "monte_carlo/include/trial_select_particle.h"

namespace feasst {

/// Select a random segment of a linear chain.
class SelectSegment : public TrialSelectParticle {
 public:
  //@{
  /** @name Arguments
   */
  
  /**
    args:
    - max_length : maximum length of selected segment. If -1 (default), then
      randomly select all possible lengths.
    - TrialSelectParticle arguments.
   */
  explicit SelectSegment(argtype args = argtype());
  explicit SelectSegment(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{
  
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
  explicit SelectSegment(std::istream& istr);
  virtual ~SelectSegment() {}

  //@}
 protected:
  void serialize_select_segment_(std::ostream& ostr) const;

 private:
  int max_length_;
};

inline std::shared_ptr<SelectSegment> MakeSelectSegment(
    argtype args = argtype()) {
  return std::make_shared<SelectSegment>(args);
}

}  // namespace feasst

#endif  // FEASST_CHAIN_SELECT_SEGMENT_H_
