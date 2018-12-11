
#ifndef FEASST_CORE_FORMULA_H_
#define FEASST_CORE_FORMULA_H_

namespace feasst {

/**
  Define a formula, used in this code for an expandable histogram.
 */
/// HWH Formula remains to be well tested, also in context of histogram.
class Formula {
 public:
  Formula& set_x0(const double x0) { x0_ = x0; return *this; }
  double x0() const { return x0_; }
  virtual double evaluate(const double x) = 0;
  virtual ~Formula() {}
 protected:
  double x0_;
};

}  // namespace feasst

#endif  // FEASST_CORE_FORMULA_H_
