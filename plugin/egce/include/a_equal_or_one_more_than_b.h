
#ifndef FEASST_EGCE_A_EQUAL_OR_ONE_MORE_THAN_B_H_
#define FEASST_EGCE_A_EQUAL_OR_ONE_MORE_THAN_B_H_

#include "flat_histogram/include/constraint.h"

namespace feasst {

/**
  Constrain the number of the first type of particle, A, to be equal to or one
  more than the number of the second type of particle, B.
 */
class AEqualOrOneMoreThanB : public Constraint {
 public:
  bool is_allowed(const System* system, const Criteria* criteria) const override {
    const int na = system->configuration().num_particles_of_type(0);
    const int nb = system->configuration().num_particles_of_type(1);
    if (na == nb || na == nb + 1) {
      return true;
    }
    return false;
  }

  virtual ~AEqualOrOneMoreThanB() {}
};

}  // namespace feasst

#endif  // FEASST_EGCE_A_EQUAL_OR_ONE_MORE_THAN_B_H_
