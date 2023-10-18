
#ifndef FEASST_STEPPERS_ANALYZE_EXAMPLE_H_
#define FEASST_STEPPERS_ANALYZE_EXAMPLE_H_

#include <vector>
#include "configuration/include/file_xyz.h"
#include "monte_carlo/include/analyze.h"

namespace feasst {

/**
  Add an analysis to FEASST by using this file as a template and instruction set.
  Follow the same 4 steps detailed in /feasst/plugin/example/include/model_example.h
  In summary, copy analyze_example.[h/cpp] to new_name.[h/cpp], replace
  AnalyzeExample with NewName, Replace ANALYZE_EXAMPLE with NEW_NAME, and finally
  cd /path/to/feasst/build; python ../py/depend-py -s ../

  For more inspiration, take a look at the other existing Analyze in the stepper
  plugin, such as Movie, Log or PairDistribution.

  In this example, compute the average geometric center of all site positions.
 */
class AnalyzeExample : public Analyze {
 public:
  //@{
  /** @name Arguments
   */

  /**
    args:
    - group_index: index of group defined in Configuration
      (default: 0, which is all existing sites).
    - Stepper arguments.

    Note that Stepper, which is the Base class of Analyze, already contains
    the vast majority of arguments that are required.
   */
  explicit AnalyzeExample(argtype args = argtype());
  explicit AnalyzeExample(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  /// Return the average geometric center of a given dimension.
  const Accumulator& geometric_center(const int dimension) const {
    return center_[dimension]; }

  /// Return the average geometric center.
  const std::vector<Accumulator>& geometric_center() const { return center_; }

  /// Write the header for the file.
  std::string header(const Criteria& criteria,
    const System& system,
    const TrialFactory& trials) const override;

  /// Size and zero the geometric center accumulator for each dimension.
  void initialize(Criteria * criteria,
      System * system,
      TrialFactory * trial_factory) override;

  /// Compute the geometric center and update the accumulator.
  void update(const Criteria& criteria,
      const System& system,
      const TrialFactory& trial_factory) override;

  /// Write the geometric center to file.
  std::string write(const Criteria& criteria,
      const System& system,
      const TrialFactory& trial_factory) override;

  // serialize
  std::string class_name() const override { return std::string("AnalyzeExample"); }
  void serialize(std::ostream& ostr) const override;
  std::shared_ptr<Analyze> create(std::istream& istr) const override {
    return std::make_shared<AnalyzeExample>(istr); }
  std::shared_ptr<Analyze> create(argtype * args) const override {
    return std::make_shared<AnalyzeExample>(args); }
  AnalyzeExample(std::istream& istr);

  //@}
 private:
  int group_index_;
  std::vector<Accumulator> center_;
};

inline std::shared_ptr<AnalyzeExample> MakeAnalyzeExample(argtype args = argtype()) {
  return std::make_shared<AnalyzeExample>(args);
}

}  // namespace feasst

#endif  // FEASST_STEPPERS_ANALYZE_EXAMPLE_H_
