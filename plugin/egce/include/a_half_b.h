
#ifndef FEASST_EGCE_A_HALF_B_H_
#define FEASST_EGCE_A_HALF_B_H_

#include "utils/include/arguments.h"
#include "monte_carlo/include/constrain_num_particles.h"

namespace feasst {

/**
  Constrain the number of the first type of particle, A, to be half of the
  number of the second type of particle, B.
 */
class AHalfB : public Constraint {
 public:
  /**
    args:
    - particle_type_A: particle type for A (default: 0).
    - particle_type_B: particle type for B (default: 1).
    - extra: allow |2*A - B| <= extra (default: 0).
   */
  AHalfB(const argtype& args = argtype());
  bool is_allowed(const System& system,
    const Criteria& criteria,
    const Acceptance& acceptance) const override;
  std::shared_ptr<Constraint> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit AHalfB(std::istream& istr);
  virtual ~AHalfB() {}

 protected:
  void serialize_a_half_b_(std::ostream& ostr) const;

 private:
  int extra_;
  ConstrainNumParticles num_A_, num_B_;
};

inline std::shared_ptr<AHalfB> MakeAHalfB(const argtype& args = argtype()) {
  return std::make_shared<AHalfB>(args);
}

}  // namespace feasst

#endif  // FEASST_EGCE_A_HALF_B_H_
