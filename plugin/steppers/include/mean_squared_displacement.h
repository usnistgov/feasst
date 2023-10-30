
#ifndef FEASST_STEPPERS_MEAN_SQUARED_DISPLACEMENT_H_
#define FEASST_STEPPERS_MEAN_SQUARED_DISPLACEMENT_H_

#include <vector>
#include "monte_carlo/include/analyze.h"

namespace feasst {

/**
  Compute the mean squared displacement of the sites over number of trials.
  Note that this is the traditional approach, not the more efficient order-n
  approach as described in Frenkel and Smit, Section 4.4.1.
  New origins are stored every fixed number of updates.
  Assume the number of particles does not change.
 */
class MeanSquaredDisplacement : public Analyze {
 public:
  //@{
  /** @name Arguments
    - updates_per_origin: set the number of updates until creation of a new
      origin (default: 1000).
    - group_index: group_index associated with configuration (default: 0).
    - Stepper arguments.
  */
  explicit MeanSquaredDisplacement(argtype args = argtype());
  explicit MeanSquaredDisplacement(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  void initialize(Criteria * criteria,
      System * system,
      TrialFactory * trial_factory) override;

  void update(const Criteria& criteria,
      const System& system,
      const TrialFactory& trial_factory) override;

  std::string write(const Criteria& criteria,
      const System& system,
      const TrialFactory& trial_factory) override;

  // serialize
  std::string class_name() const override {
    return std::string("MeanSquaredDisplacement"); }
  void serialize(std::ostream& ostr) const override;
  std::shared_ptr<Analyze> create(std::istream& istr) const override {
    return std::make_shared<MeanSquaredDisplacement>(istr); }
  std::shared_ptr<Analyze> create(argtype * args) const override {
    return std::make_shared<MeanSquaredDisplacement>(args); }
  explicit MeanSquaredDisplacement(std::istream& istr);

  //@}
 private:
  int updates_since_origin_;
  int updates_per_origin_;
  int group_index_;
  std::vector<Select> origins_;
  std::vector<Accumulator> msd_;

  int num_frames_() const { return static_cast<int>(msd_.size()); }

  void update_msd_(const int updates,
      const Configuration& config,
      const Select& old_parts);
};

inline std::shared_ptr<MeanSquaredDisplacement> MakeMeanSquaredDisplacement(
    argtype args = argtype()) {
  return std::make_shared<MeanSquaredDisplacement>(args);
}

}  // namespace feasst

#endif  // FEASST_STEPPERS_MEAN_SQUARED_DISPLACEMENT_H_
