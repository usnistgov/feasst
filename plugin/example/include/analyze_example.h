
#ifndef FEASST_STEPPERS_ANALYZE_EXAMPLE_H_
#define FEASST_STEPPERS_ANALYZE_EXAMPLE_H_

#include <memory>
#include <map>
#include <string>
#include <vector>
#include "monte_carlo/include/analyze.h"

namespace feasst {

class Accumulator;
class AveragePosition;
class MonteCarlo;
class VisitConfiguration;

typedef std::map<std::string, std::string> argtype;

/**
  Add an analysis to FEASST by using this file as a template and instruction set.
  Follow the same 4 steps detailed in /feasst/plugin/example/include/model_example.h
  In summary, copy analyze_example.[h/cpp] to new_name.[h/cpp], replace
  AnalyzeExample with NewName, Replace ANALYZE_EXAMPLE with NEW_NAME, and finally
  cd /path/to/feasst/build; python ../dev/tools/depend-py -s ../

  For more inspiration, take a look at the other existing Analyze in the stepper
  plugin, such as Movie, Log or PairDistribution.

  In this example, compute the average geometric center of all site positions.
 */
class AnalyzeExample : public Analyze {
 public:
  //@{
  /** @name Arguments
    - group: name of group defined in Configuration (default: all sites).
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
  const Accumulator& geometric_center(const int dimension) const;

  /// Return the average geometric center.
  const std::vector<std::unique_ptr<Accumulator> >& geometric_center() const;

  /// Write the header for the file.
  std::string header(const MonteCarlo& mc) const override;

  /// Size and zero the geometric center accumulator for each dimension.
  void initialize(MonteCarlo * mc) override;

  /// Compute the geometric center and update the accumulator.
  void update(const MonteCarlo& mc) override;

  /// Write the geometric center to file.
  std::string write(const MonteCarlo& mc) override;

  // serialize
  std::string class_name() const override {
    return std::string("AnalyzeExample"); }
  void serialize(std::ostream& ostr) const override;
  std::shared_ptr<Analyze> create(std::istream& istr) const override {
    return std::make_shared<AnalyzeExample>(istr); }
  std::shared_ptr<Analyze> create(argtype * args) const override {
    return std::make_shared<AnalyzeExample>(args); }
  explicit AnalyzeExample(std::istream& istr);
  virtual ~AnalyzeExample();

  //@}
 private:
  int group_index_ = 0;
  std::vector<std::unique_ptr<Accumulator> > center_;
  std::unique_ptr<AveragePosition> loop_;
  std::unique_ptr<VisitConfiguration> visit_;

  // not serialized, assumes initialize is run atleast once
  std::string group_;
};

}  // namespace feasst

#endif  // FEASST_STEPPERS_ANALYZE_EXAMPLE_H_
