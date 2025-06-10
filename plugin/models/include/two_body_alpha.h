
#ifndef FEASST_MODELS_TWO_BODY_ALPHA_H_
#define FEASST_MODELS_TWO_BODY_ALPHA_H_

#include "configuration/include/model_params.h"
#include "system/include/model_two_body.h"

namespace feasst {

typedef std::map<std::string, std::string> argtype;

/**
  The two-body alpha model is defined as

  \f$U = s\epsilon \left(\frac{\sigma}{r}\right)^\alpha\f$

  where \f$s\f$ is the prefactor and \f$\alpha\f$ is the power,
  both of which are given as class-specific parameters.
  Note that \f$s\f$ is required for negative values because \f$\epsilon\f$
  cannot be negative due to the standard mixing rules.
 */
class TwoBodyAlpha : public ModelTwoBody {
 public:
  /**
    args:
    - alpha: set the value of \f$\alpha\f$.
      A comma-separated list provides an alpha for each potential.
    - s: set the value of \f$s\f$ with a comma-separated list.
   */
  explicit TwoBodyAlpha(argtype args = argtype());
  explicit TwoBodyAlpha(argtype * args);

  double energy(
    const double squared_distance,
    const int type1,
    const int type2,
    const ModelParams& model_params) override;

  std::shared_ptr<Model> create(std::istream& istr) const override {
    return std::make_shared<TwoBodyAlpha>(istr); }
  std::shared_ptr<Model> create(argtype * args) const override {
    return std::make_shared<TwoBodyAlpha>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit TwoBodyAlpha(std::istream& istr);
  virtual ~TwoBodyAlpha() {}

 protected:
  void serialize_two_body_alpha_(std::ostream& ostr) const;

 private:
  std::vector<double> alpha_;
  std::vector<double> s_;
};

inline std::shared_ptr<TwoBodyAlpha> MakeTwoBodyAlpha(
    argtype args = argtype()) {
  return std::make_shared<TwoBodyAlpha>(args);
}

}  // namespace feasst

#endif  // FEASST_MODELS_TWO_BODY_ALPHA_H_
