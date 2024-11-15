
#ifndef FEASST_STEPPERS_WRAP_PARTICLES_H_
#define FEASST_STEPPERS_WRAP_PARTICLES_H_

#include "monte_carlo/include/modify.h"

namespace feasst {

/**
  Wrap particles into the periodic domains based on the position of the first
  site.
 */
class WrapParticles : public ModifyUpdateOnly {
 public:
  //@{
  /** @name Arguments
    - Stepper arguments.
   */
  explicit WrapParticles(argtype args = argtype());
  explicit WrapParticles(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  void update(MonteCarlo * mc) override;

  std::string class_name() const override { return std::string("WrapParticles"); }

  // serialize
  void serialize(std::ostream& ostr) const override;
  std::shared_ptr<Modify> create(std::istream& istr) const override {
    return std::make_shared<WrapParticles>(istr); }
  std::shared_ptr<Modify> create(argtype * args) const override {
    return std::make_shared<WrapParticles>(args); }
  explicit WrapParticles(std::istream& istr);
  //@}
};

inline std::shared_ptr<WrapParticles> MakeWrapParticles(
    argtype args = argtype()) {
  return std::make_shared<WrapParticles>(args);
}

}  // namespace feasst

#endif  // FEASST_STEPPERS_WRAP_PARTICLES_H_
