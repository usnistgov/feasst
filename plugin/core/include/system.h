
#ifndef FEASST_CORE_SYSTEM_H_
#define FEASST_CORE_SYSTEM_H_

#include <vector>
#include "core/include/debug.h"
#include "core/include/configuration.h"
#include "core/include/visit_model.h"
#include "core/include/model_lj.h"

namespace feasst {

class System {
 public:
  int dimension() const { return configurations_[0].dimension(); }

  void add_configuration(Configuration configuration) { configurations_.push_back(configuration); }
  Configuration* configuration(const int iConfiguration) { return &configurations_[0]; }
  const Configuration& config() const { return configurations_[0]; }
  int num_configurations() const { return configurations_.size(); }

//  void add(OneBody one_body) { one_bodies_.push_back(one_body); }
//  void add(TwoBody two_body) { two_bodies_.push_back(two_body); }
//  TwoBody* pair(const int iTwoBody) { return &two_bodies_[iTwoBody]; }
//  unsigned int num_two_bodies() const { return two_bodies_.size(); }

  int num_particles() const {
    int num = 0;
    for (int iConf = 0; iConf < num_configurations(); ++iConf) {
      num += configurations_[iConf].num_particles();
    }
    return num;
  }

  Configuration * configurationByPart(const int iPart) {
    ASSERT(configurations_.size() == 1, "error");
    return &configurations_[0];
  }

  Particle particle(const int particle_index) {
    return configurationByPart(particle_index)->particle(particle_index);
  }

  void default_system() {
    Configuration config;
    config.default_configuration();
    add_configuration(config);
  }

  double energy() {
    return energy(model_);
  }

  double energy(const Select& selection) {
    return model_.compute(visit_model_, configurations_[0], selection);
  }

  double energy_of_selection() {
    return model_.compute(visit_model_, configurations_[0], configurations_[0].selection_of_all());
  }

  double energy(const ModelTwoBody& model) {
    ASSERT(configurations_.size() > 0, "size error");
    visit_model_.compute(configurations_[0], model);
    return visit_model_.energy();
  }

  ModelTwoBody* model() { return &model_; }
  VisitModel* visitor() { return &visit_model_; }

 private:
  std::vector<Configuration> configurations_;
  VisitModel visit_model_;
//  vector<TwoBody> two_bodies_;
//  vector<OneBody> one_bodies_;
  ModelLJ model_;
//  vector<ModelTwoBody> two_body_models_;
//  vector<ModelOneBody> one_body_models_;
};

}  // namespace feasst

#endif  // FEASST_CORE_SYSTEM_H_
