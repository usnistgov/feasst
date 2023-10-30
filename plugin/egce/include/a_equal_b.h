
#ifndef FEASST_EGCE_A_EQUAL_B_H_
#define FEASST_EGCE_A_EQUAL_B_H_

#include "utils/include/arguments.h"
#include "monte_carlo/include/constrain_num_particles.h"

namespace feasst {

/**
  Constrain the number of the first type of particle, A, to be equal to the
  number of the second type of particle, B.
 */
class AEqualB : public Constraint {
 public:
  //@{
  /** @name Arguments
    - particle_type_A: particle type for A (default: 0).
    - particle_type_B: particle type for B (default: 1).
    - extra_A: allow this many extra A particles (default: 0).
   */
  explicit AEqualB(argtype args = argtype());
  explicit AEqualB(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  bool is_allowed(const System& system,
    const Criteria& criteria,
    const Acceptance& acceptance) const override;
  std::shared_ptr<Constraint> create(std::istream& istr) const override {
    return std::make_shared<AEqualB>(istr); }
  std::shared_ptr<Constraint> create(argtype * args) const override {
    return std::make_shared<AEqualB>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit AEqualB(std::istream& istr);
  virtual ~AEqualB() {}

  //@}
 protected:
  void serialize_a_equal_b_(std::ostream& ostr) const;

 private:
  int extra_A_;
  ConstrainNumParticles num_A_, num_B_;
};

inline std::shared_ptr<AEqualB> MakeAEqualB(argtype args = argtype()) {
  return std::make_shared<AEqualB>(args);
}

}  // namespace feasst

#endif  // FEASST_EGCE_A_EQUAL_B_H_
