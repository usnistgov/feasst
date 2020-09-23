
#ifndef FEASST_STEPPERS_VOLUME_H_
#define FEASST_STEPPERS_VOLUME_H_

#include "monte_carlo/include/analyze.h"
#include "math/include/accumulator.h"

namespace feasst {

/**
  Accumulate average volume.
 */
class Volume : public Analyze {
 public:
  explicit Volume(const argtype &args = argtype());

  std::string header(const Criteria& criteria,
    const System& system,
    const TrialFactory& trials) const override;

  void initialize(Criteria * criteria,
      System * system,
      TrialFactory * trial_factory) override;

  void update(const Criteria& criteria,
      const System& system,
      const TrialFactory& trial_factory) override;

  std::string write(const Criteria& criteria,
      const System& system,
      const TrialFactory& trial_factory) override;

  const Accumulator& energy() const { return accumulator(); }

  // serialize
  std::string class_name() const override { return std::string("Volume"); }
  std::shared_ptr<Analyze> create(std::istream& istr) const override {
    return std::make_shared<Volume>(istr); }
  void serialize(std::ostream& ostr) const override;
  explicit Volume(std::istream& istr);
};

inline std::shared_ptr<Volume> MakeVolume(const argtype &args = argtype()) {
  return std::make_shared<Volume>(args);
}

}  // namespace feasst

#endif  // FEASST_STEPPERS_VOLUME_H_
