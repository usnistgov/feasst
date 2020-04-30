
#ifndef FEASST_EGCE_A_TWICE_B_H_
#define FEASST_EGCE_A_TWICE_B_H_

#include "monte_carlo/include/constraint.h"

namespace feasst {

/**
  Constrain the number of the first type of particle, A, to be twice that of
  the number of the second type of particle, B.
  But also allow a slight deviation, |na/2 - nb| <= 0.5
 */
class ATwiceB : public Constraint {
 public:
  ATwiceB() { class_name_ = "ATwiceB"; }
  bool is_allowed(const System* system,
    const Criteria* criteria,
    const Acceptance& acceptance) const override;
  std::shared_ptr<Constraint> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit ATwiceB(std::istream& istr);
  virtual ~ATwiceB() {}

 protected:
  void serialize_a_twice_b_(std::ostream& ostr) const;
};

inline std::shared_ptr<ATwiceB> MakeATwiceB() {
  return std::make_shared<ATwiceB>();
}

}  // namespace feasst

#endif  // FEASST_EGCE_A_TWICE_B_H_
