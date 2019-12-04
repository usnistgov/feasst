
#ifndef FEASST_STEPPERS_WALL_CLOCK_LIMIT_H_
#define FEASST_STEPPERS_WALL_CLOCK_LIMIT_H_

#include "monte_carlo/include/modify.h"

namespace feasst {

/**
  Terminate the simulation after a given number of CPU hours in order to
  prevent fragmented checkpoint files.
 */
class WallClockLimit : public ModifyUpdateOnly {
 public:
  /**
    args:
    - max_hours: maximum number of wall clock hours until job termination.
   */
  WallClockLimit(const argtype &args = argtype());

  void update(Criteria * criteria,
      System * system,
      TrialFactory * trial_factory) override;

  // serialize
  std::string class_name() const override { return std::string("WallClockLimit"); }
  void serialize(std::ostream& ostr) const override;
  std::shared_ptr<Modify> create(std::istream& istr) const override {
    return std::make_shared<WallClockLimit>(istr); }
  WallClockLimit(std::istream& istr);

 private:
  double max_hours_ = 0;
};

inline std::shared_ptr<WallClockLimit> MakeWallClockLimit(
    const argtype &args = argtype()) {
  return std::make_shared<WallClockLimit>(args);
}

}  // namespace feasst

#endif  // FEASST_STEPPERS_WALL_CLOCK_LIMIT_H_
