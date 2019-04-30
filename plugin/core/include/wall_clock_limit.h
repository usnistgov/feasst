
#ifndef FEASST_CORE_WALL_CLOCK_LIMIT_H_
#define FEASST_CORE_WALL_CLOCK_LIMIT_H_

#include "core/include/modify.h"

namespace feasst {

/**
 */
class WallClockLimit : public ModifyUpdateOnly {
 public:
  /**
    max_hours : maximum number of wall clock hours until job termination.
   */
  WallClockLimit(const argtype &args = argtype()) : ModifyUpdateOnly(args) {
    args_.init(args);
    max_hours_ = args_.key("max_hours").dble();
    set_steps_per(args_.key("steps_per").dflt("1").integer());
  }
  void update(std::shared_ptr<Criteria> criteria,
      System * system,
      TrialFactory * trial_factory) override {
    const double hours = double(clock())/double(CLOCKS_PER_SEC)/60./60.;
    ASSERT(hours < max_hours_, "wall clock hours(" << hours << ") exceed " <<
      "the maximum(" << max_hours_ << ")");
  }

  void serialize(std::ostream& ostr) const override {
    ostr << class_name_ << " ";
    feasst_serialize_version(1, ostr);
    feasst_serialize(max_hours_, ostr);
  }

  std::shared_ptr<Modify> create(std::istream& istr) const override {
    feasst_deserialize_version(istr);
    auto modify = std::make_shared<WallClockLimit>();
    feasst_deserialize(&(modify->max_hours_), istr);
    return modify;
  }

 private:
  const std::string class_name_ = "WallClockLimit";
  double max_hours_ = 0;
};

inline std::shared_ptr<WallClockLimit> MakeWallClockLimit(const argtype &args = argtype()) {
  return std::make_shared<WallClockLimit>(args);
}

}  // namespace feasst

#endif  // FEASST_CORE_WALL_CLOCK_LIMIT_H_
