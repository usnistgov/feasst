
#ifndef FEASST_SYSTEM_MODEL_ONE_BODY_H_
#define FEASST_SYSTEM_MODEL_ONE_BODY_H_

#include "system/include/model.h"
#include "system/include/visit_model.h"

namespace feasst {

class Site;
class Configuration;

class ModelOneBody : public Model {
 public:
  double compute(
      const ModelParams& model_params,
      const int group_index,
      Configuration * config,
      VisitModel * visitor) const override {
    visitor->compute(*this, model_params, config, group_index);
    return visitor->energy();
  }
  double compute(
      const int group_index,
      Configuration * config,
      VisitModel * visitor) const override {
    visitor->compute(*this, config, group_index);
    return visitor->energy();
  }
  double compute(
      const ModelParams& model_params,
      Configuration * config,
      VisitModel * visitor) const override {
    visitor->compute(*this, model_params, config, 0);
    return visitor->energy();
  }
  double compute(
      Configuration * config,
      VisitModel * visitor) const override {
    visitor->compute(*this, config, 0);
    return visitor->energy();
  }
  double compute(
      const ModelParams& model_params,
      const Select& selection,
      const int group_index,
      Configuration * config,
      VisitModel * visitor) const override {
    visitor->compute(*this, model_params, selection, config, group_index);
    return visitor->energy();
  }
  double compute(
      const Select& selection,
      const int group_index,
      Configuration * config,
      VisitModel * visitor) const override {
    visitor->compute(*this, selection, config, group_index);
    return visitor->energy();
  }
  double compute(
      const ModelParams& model_params,
      const Select& selection,
      Configuration * config,
      VisitModel * visitor) const override {
    visitor->compute(*this, model_params, selection, config, 0);
    return visitor->energy();
  }
  double compute(
      const Select& selection,
      Configuration * config,
      VisitModel * visitor) const override {
    visitor->compute(*this, selection, config, 0);
    return visitor->energy();
  }

  /// Return the energy given the wrapped coordinates, site, config and params.
  virtual double energy(
    const Position& wrapped_site,
    const Site& site,
    const Configuration& config,
    const ModelParams& model_params) const = 0;

  /// Same as above, but assume that the site position is already wrapped.
  double energy_no_wrap(
    const Site& site,
    const Configuration& config,
    const ModelParams& model_params) const;

  virtual ~ModelOneBody() {}

 private:
  using Model::energy;  // remove hidden overloaded virtual function warnings.
};

}  // namespace feasst

#endif  // FEASST_SYSTEM_MODEL_ONE_BODY_H_
