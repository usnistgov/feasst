
#ifndef FEASST_SYSTEM_MODEL_ONE_BODY_H_
#define FEASST_SYSTEM_MODEL_ONE_BODY_H_

#include "system/include/model.h"

namespace feasst {

class Configuration;
class Position;
class Site;
class VisitModel;

class ModelOneBody : public Model {
 public:
  ModelOneBody() { class_name_ = "ModelOneBody"; }

  double compute(
    const ModelParams& model_params,
    const int group_index,
    Configuration * config,
    VisitModel * visitor) override;
  double compute(
    const int group_index,
    Configuration * config,
    VisitModel * visitor) override;
  double compute(
    const ModelParams& model_params,
    Configuration * config,
    VisitModel * visitor) override;
  double compute(
    Configuration * config,
    VisitModel * visitor) override;
  double compute(
    const ModelParams& model_params,
    const Select& selection,
    const int group_index,
    Configuration * config,
    VisitModel * visitor) override;
  double compute(
    const Select& selection,
    const int group_index,
    Configuration * config,
    VisitModel * visitor) override;
  double compute(
    const ModelParams& model_params,
    const Select& selection,
    Configuration * config,
    VisitModel * visitor) override;
  double compute(
    const Select& selection,
    Configuration * config,
    VisitModel * visitor) override;

  int num_body() const override { return 1; }

  /// Return the energy given the wrapped coordinates, site, config and params.
  virtual double energy(
    const Position& wrapped_site,
    const Site& site,
    const Configuration& config,
    const ModelParams& model_params) = 0;

  /// Same as above, but assume that the site position is already wrapped.
  double energy_no_wrap(
    const Site& site,
    const Configuration& config,
    const ModelParams& model_params);

  virtual ~ModelOneBody() {}
  explicit ModelOneBody(std::istream& istr) : Model(istr) {}

 private:
  using Model::energy;  // remove hidden overloaded virtual function warnings.
};

}  // namespace feasst

#endif  // FEASST_SYSTEM_MODEL_ONE_BODY_H_
