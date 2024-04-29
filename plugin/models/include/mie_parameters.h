
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
  The mie_ideal_deviation parameter has the default mixing rule:
  \f$ \lambda_{ij} = 0 \f$
  and must be input manually.
 */
class MieIdealDeviation : public ModelParam {
 public:
  MieIdealDeviation() : ModelParam() { class_name_ = "mie_ideal_deviation"; }
  std::shared_ptr<ModelParam> create(std::istream& istr) const override {
    return std::make_shared<MieIdealDeviation>(istr); }
  void serialize(std::ostream& ostr) const override;
  explicit MieIdealDeviation(std::istream& istr);
  virtual ~MieIdealDeviation() {}
 private:
  double mix_(const double value1, const double value2) override;
};

inline std::shared_ptr<MieIdealDeviation> MakeMieIdealDeviation() {
  return std::make_shared<MieIdealDeviation>();
}

}  // namespace feasst

#endif  // FEASST_MODELS_MIE_PARAMETERS_H_
