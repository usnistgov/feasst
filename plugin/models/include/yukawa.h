
#ifndef FEASST_MODELS_YUKAWA_H_
#define FEASST_MODELS_YUKAWA_H_

#include <sstream>
#include "system/include/model_two_body.h"

namespace feasst {

/**
  The screened Coulomb (Yukawa) potential is given by
  \f$ U = \epsilon \exp(-\kappa (r/\sigma - 1)) / (r/\sigma) \f$
  Note that \f$\kappa\f$ is dimensionless in this definition.
 */
class Yukawa : public ModelTwoBody {
 public:
  //@{
  /** @name Arguments
    - kappa: set the value of the kappa parameter (default: 1).
   */
  explicit Yukawa(argtype args = argtype());
  explicit Yukawa(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  double energy(
      const double squared_distance,
      const int type1,
      const int type2,
      const ModelParams& model_params) override;

  /// Set the value of the kappa parameter.
  void set_kappa(const double kappa = 1) { kappa_ = kappa; }
  double kappa() const { return kappa_; }

  std::shared_ptr<Model> create(std::istream& istr) const override {
    return std::make_shared<Yukawa>(istr); }
  std::shared_ptr<Model> create(argtype * args) const override {
    return std::make_shared<Yukawa>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit Yukawa(std::istream& istr);
  virtual ~Yukawa() {}

  //@}
 private:
  double kappa_;
};

inline std::shared_ptr<Yukawa> MakeYukawa(argtype args = argtype()) {
  return std::make_shared<Yukawa>(args);
}

}  // namespace feasst

#endif  // FEASST_MODELS_YUKAWA_H_
