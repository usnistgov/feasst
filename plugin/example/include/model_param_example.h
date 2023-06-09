
#ifndef FEASST_EXAMPLE_MODEL_PARAM_EXAMPLE_H_
#define FEASST_EXAMPLE_MODEL_PARAM_EXAMPLE_H_

#include <memory>
#include "configuration/include/model_params.h"

namespace feasst {

/**
  Add a ModelParam to FEASST by using this file as a template and instruction set.
  Follow the same steps detailed in /feasst/plugin/example/README.rst.
  In summary, copy model_param_example.[h/cpp] to new_name.[h/cpp], replace
  Gamma with NewName, then replace MODEL_PARAM_EXAMPLE with NEW_NAME.

  While sigma, epsilon and cutoff are already implemented in Model,
  and the lambda parameter was already implemented in LennardJonesAlpha,
  the gamma parameter is not yet defined.
  Thus, we will go through the exercise here.
 */
class Gamma : public ModelParam {
 public:
  Gamma() : ModelParam() { class_name_ = "gamma"; }
  std::shared_ptr<ModelParam> create(std::istream& istr) const override {
    return std::make_shared<Gamma>(istr); }
  void serialize(std::ostream& ostr) const override;
  explicit Gamma(std::istream& istr);
  virtual ~Gamma() {}
 private:
  double mix_(const double value1, const double value2) override;
};

}  // namespace feasst

#endif  // FEASST_EXAMPLE_MODEL_PARAM_EXAMPLE_H_
