
#ifndef FEASST_MODELS_MIE_PARAMETERS_H_
#define FEASST_MODELS_MIE_PARAMETERS_H_

#include "configuration/include/model_params.h"

namespace feasst {

/**
  The mie_lambda_r parameter has the default mixing rule:
  \f$ \lambda_{r,ij} - 3 = \sqrt{(\lambda_{r,ii} - 3)(\lambda_{r,jj}-3)} \f$
 */
class MieLambdaR : public ModelParam {
 public:
  MieLambdaR() : ModelParam() { class_name_ = "mie_lambda_r"; }
  std::shared_ptr<ModelParam> create(std::istream& istr) const override {
    return std::make_shared<MieLambdaR>(istr); }
  void serialize(std::ostream& ostr) const override;
  explicit MieLambdaR(std::istream& istr);
  virtual ~MieLambdaR() {}
 private:
  double mix_(const double value1, const double value2) override;
};

inline std::shared_ptr<MieLambdaR> MakeMieLambdaR() {
  return std::make_shared<MieLambdaR>();
}

/**
  The mie_lambda_a parameter has the default mixing rule:
  \f$ \lambda_{a,ij} - 3 = \sqrt{(\lambda_{a,ii} - 3)(\lambda_{a,jj}-3)} \f$
 */
class MieLambdaA : public ModelParam {
 public:
  MieLambdaA() : ModelParam() { class_name_ = "mie_lambda_a"; }
  std::shared_ptr<ModelParam> create(std::istream& istr) const override {
    return std::make_shared<MieLambdaA>(istr); }
  void serialize(std::ostream& ostr) const override;
  explicit MieLambdaA(std::istream& istr);
  virtual ~MieLambdaA() {}
 private:
  double mix_(const double value1, const double value2) override;
};

inline std::shared_ptr<MieLambdaA> MakeMieLambdaA() {
  return std::make_shared<MieLambdaA>();
}

/**
  The mie_prefactor parameter is given by

  \f$ \epsilon n/(n-m)(n/m)^{m/(n-m} \f$

  where n and m are mie_lambda_r and mie_lambda_a, respectively.
  These are precomputed and do not need to be entered manually but are instead
  precomputed and stored for optimization.
 */
class MiePrefactor : public ModelParam {
 public:
  MiePrefactor() : ModelParam() { class_name_ = "mie_prefactor"; }
  double compute(const int type1, const int type2,
    const ModelParams& model_params) override;
  std::shared_ptr<ModelParam> create(std::istream& istr) const override {
    return std::make_shared<MiePrefactor>(istr); }
  void serialize(std::ostream& ostr) const override;
  explicit MiePrefactor(std::istream& istr);
  virtual ~MiePrefactor() {}
};

///**
//  The mie_ideal_deviation parameter has the default mixing rule:
//  \f$ \lambda_{ij} = 0 \f$
//  and must be input manually.
// */
//class MieIdealDeviation : public ModelParam {
// public:
//  MieIdealDeviation() : ModelParam() { class_name_ = "mie_ideal_deviation"; }
//  std::shared_ptr<ModelParam> create(std::istream& istr) const override {
//    return std::make_shared<MieIdealDeviation>(istr); }
//  void serialize(std::ostream& ostr) const override;
//  explicit MieIdealDeviation(std::istream& istr);
//  virtual ~MieIdealDeviation() {}
// private:
//  double mix_(const double value1, const double value2) override;
//};
//
//inline std::shared_ptr<MieIdealDeviation> MakeMieIdealDeviation() {
//  return std::make_shared<MieIdealDeviation>();
//}

}  // namespace feasst

#endif  // FEASST_MODELS_MIE_PARAMETERS_H_
