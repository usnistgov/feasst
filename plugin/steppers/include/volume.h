
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
  //@{
  /** @name Arguments
    - Stepper arguments.
   */
  explicit Volume(argtype args = argtype());
  explicit Volume(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  std::string header(const MonteCarlo& mc) const override;
  void initialize(MonteCarlo * mc) override;
  void update(const MonteCarlo& mc) override;
  std::string write(const MonteCarlo& mc) override;

  const Accumulator& volume() const { return accumulator(); }

  // serialize
  std::string class_name() const override { return std::string("Volume"); }
  std::shared_ptr<Analyze> create(std::istream& istr) const override {
    return std::make_shared<Volume>(istr); }
  std::shared_ptr<Analyze> create(argtype * args) const override {
    return std::make_shared<Volume>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit Volume(std::istream& istr);
  //@}
};

inline std::shared_ptr<Volume> MakeVolume(argtype args = argtype()) {
  return std::make_shared<Volume>(args);
}

}  // namespace feasst

#endif  // FEASST_STEPPERS_VOLUME_H_
