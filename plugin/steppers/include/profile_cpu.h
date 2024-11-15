
#ifndef FEASST_STEPPERS_PROFILE_CPU_H_
#define FEASST_STEPPERS_PROFILE_CPU_H_

#include "monte_carlo/include/analyze.h"

namespace feasst {

/**
  Periodically write the profile of where CPU time is spent.
 */
class ProfileCPU : public AnalyzeWriteOnly {
 public:
  //@{
  /** @name Arguments
    - Stepper arguments.
   */
  explicit ProfileCPU(argtype args = argtype());
  explicit ProfileCPU(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  std::string header(const MonteCarlo& mc) const override;
  void initialize(MonteCarlo * mc) override;
  std::string write(const MonteCarlo& mc) override;

  // serialize
  std::string class_name() const override { return std::string("ProfileCPU"); }
  std::shared_ptr<Analyze> create(std::istream& istr) const override {
    return std::make_shared<ProfileCPU>(istr); }
  std::shared_ptr<Analyze> create(argtype * args) const override {
    return std::make_shared<ProfileCPU>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit ProfileCPU(std::istream& istr);
  //@}
};

inline std::shared_ptr<ProfileCPU> MakeProfileCPU(argtype args = argtype()) {
  return std::make_shared<ProfileCPU>(args);
}

}  // namespace feasst

#endif  // FEASST_STEPPERS_PROFILE_CPU_H_
