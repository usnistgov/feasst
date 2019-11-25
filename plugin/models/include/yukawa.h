
#ifndef FEASST_MODELS_YUKAWA_H_
#define FEASST_MODELS_YUKAWA_H_

#include "system/include/model_two_body.h"
#include "math/include/constants.h"

namespace feasst {

/**
  The screened Coulomb (Yukawa) potential is given by
  \f$ U = \epsilon \exp(-\kappa (r/\sigma - 1)) / (r/\sigma) \f$
  Note that \f$\kappa\f$ is dimensionless in this definition.
 */
class Yukawa : public ModelTwoBody {
 public:
  Yukawa();

  double energy(
      const double squared_distance,
      const int type1,
      const int type2,
      const ModelParams& model_params) const override {
    const double epsilon = model_params.mixed_epsilon()[type1][type2];
    const double sigma = model_params.mixed_sigma()[type1][type2];
    const double distance = sqrt(squared_distance);
    TRACE("epsilon " << epsilon << " distance " << distance << " kappa " << kappa_);
    return epsilon*exp(-kappa_*(distance/sigma - 1.))/(distance/sigma);
  }

  /// Set the value of the kappa parameter.
  void set_kappa(const double kappa = 1) { kappa_ = kappa; }
  double kappa() const { return kappa_; }

  std::shared_ptr<Model> create(std::istream& istr) const override {
    auto model = std::make_shared<Yukawa>();
    int version;
    istr >> version;
    double kappa;
    istr >> kappa;
    model->set_kappa(kappa);
    return model;
  }

  void serialize(std::ostream& ostr) const override {
    ostr << class_name_ << " "
         << "1 " // version
         << kappa_ << " ";
  }

  virtual ~Yukawa() {}

 private:
  const std::string class_name_ = "Yukawa";
  double kappa_;
};

inline std::shared_ptr<Yukawa> MakeYukawa() {
  return std::make_shared<Yukawa>();
}

}  // namespace feasst

#endif  // FEASST_MODELS_YUKAWA_H_
