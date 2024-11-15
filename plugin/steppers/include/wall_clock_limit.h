
#ifndef FEASST_STEPPERS_WALL_CLOCK_LIMIT_H_
#define FEASST_STEPPERS_WALL_CLOCK_LIMIT_H_

#include "monte_carlo/include/analyze.h"

namespace feasst {

/**
  Terminate the simulation after a given number of CPU hours in order to
  prevent fragmented checkpoint files.
  Thanks to Joshua Anderson for suggesting this.
 */
class WallClockLimit : public AnalyzeUpdateOnly {
 public:
  //@{
  /** @name Arguments
    - max_hours: maximum number of wall clock hours until job termination.
    - Stepper arguments.
   */
  explicit WallClockLimit(argtype args = argtype());
  explicit WallClockLimit(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  void update(const MonteCarlo& mc) override;

  // serialize
  std::string class_name() const override { return std::string("WallClockLimit"); }
  void serialize(std::ostream& ostr) const override;
  std::shared_ptr<Analyze> create(std::istream& istr) const override {
    return std::make_shared<WallClockLimit>(istr); }
  std::shared_ptr<Analyze> create(argtype * args) const override {
    return std::make_shared<WallClockLimit>(args); }
  explicit WallClockLimit(std::istream& istr);

  //@}
 private:
  double max_hours_ = 0;
};

inline std::shared_ptr<WallClockLimit> MakeWallClockLimit(
    argtype args = argtype()) {
  return std::make_shared<WallClockLimit>(args);
}

}  // namespace feasst

#endif  // FEASST_STEPPERS_WALL_CLOCK_LIMIT_H_
