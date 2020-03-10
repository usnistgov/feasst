
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
  AEqualOrOneMoreThanB() { class_name_ = "AEqualOrOneMoreThanB"; }
  bool is_allowed(const System* system,
    const Criteria* criteria,
    const Acceptance& acceptance) const override;
  std::shared_ptr<Constraint> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit AEqualOrOneMoreThanB(std::istream& istr);
  virtual ~AEqualOrOneMoreThanB() {}

 protected:
  void serialize_a_equal_or_one_more_than_b_(std::ostream& ostr) const;
};

inline std::shared_ptr<AEqualOrOneMoreThanB> MakeAEqualOrOneMoreThanB() {
  return std::make_shared<AEqualOrOneMoreThanB>();
}

}  // namespace feasst

#endif  // FEASST_EGCE_A_EQUAL_OR_ONE_MORE_THAN_B_H_
