
#ifndef FEASST_MODELS_SQUARE_WELL_H_
#define FEASST_MODELS_SQUARE_WELL_H_

#include "system/include/model_two_body.h"

namespace feasst {

/**
  The square well model is defined as

  \f$U(r) = \left\{
    \begin{array}{lr}
      \infty & r < \sigma \\
      -\epsilon & \sigma \le r \le r_c \\
      0 & r > r_c
    \end{array}
  \right.\f$
 */
class SquareWell : public ModelTwoBody {
 public:
  SquareWell() { class_name_ = "SquareWell"; }

  double energy(
    const double squared_distance,
    const int type1,
    const int type2,
    const ModelParams& model_params) override;

  std::shared_ptr<Model> create(std::istream& istr) const override {
    return std::make_shared<SquareWell>(istr); }
  std::shared_ptr<Model> create(argtype * args) const override {
    return std::make_shared<SquareWell>(); }
  void serialize(std::ostream& ostr) const override;
  explicit SquareWell(std::istream& istr);
  virtual ~SquareWell() {}
};

inline std::shared_ptr<SquareWell> MakeSquareWell() {
  return std::make_shared<SquareWell>();
}

}  // namespace feasst

#endif  // FEASST_MODELS_SQUARE_WELL_H_
