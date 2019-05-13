
#ifndef FEASST_SYSTEM_MODEL_TWO_BODY_H_
#define FEASST_SYSTEM_MODEL_TWO_BODY_H_

#include "system/include/model.h"
#include "system/include/visit_model.h"
#include "configuration/include/particle.h"
#include "math/include/constants.h"

namespace feasst {

class ModelTwoBody : public Model {
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
  virtual double energy(
      const double squared_distance,
      const int type1,
      const int type2,
      const ModelParams& model_params) const = 0;
  virtual ~ModelTwoBody() {}
};

class ModelTwoBodyFactory : public ModelTwoBody {
 public:
  void add_model(std::shared_ptr<ModelTwoBody> model) {
    models_.push_back(model); }
  const std::vector<std::shared_ptr<ModelTwoBody> >& models() const {
    return models_; }

  double energy(
      const double squared_distance,
      const int type1,
      const int type2,
      const ModelParams& model_params) const override {
    double en = 0;
    for (const std::shared_ptr<ModelTwoBody> model : models_) {
      en += model->energy(squared_distance, type1, type2, model_params);
    }
//    int index = 0;
//    while ((index < static_cast<int>(models_.size())) &&
//           (en < NEAR_INFINITY)) {
//      en += models_[index]->energy(squared_distance,
//                                   type1, type2, model_params);
//      ++index;
//    }
    //INFO("en " << en);
    return en;
  }

  void precompute(const ModelParams& existing) override {
    for (const std::shared_ptr<ModelTwoBody> model : models_) {
      model->precompute(existing);
    }
  }

 private:
  std::vector<std::shared_ptr<ModelTwoBody> > models_;
};

}  // namespace feasst

#endif  // FEASST_SYSTEM_MODEL_TWO_BODY_H_
