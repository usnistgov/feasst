
#ifndef FEASST_CORE_VISIT_MODEL_H_
#define FEASST_CORE_VISIT_MODEL_H_

#include "core/include/model.h"
#include "core/include/configuration.h"

namespace feasst {

class ModelOneBody;
class ModelTwoBody;

class VisitModel {
 public:
  // hwh depreciate
  virtual void loop_by_particle(const Configuration& config, const ModelOneBody& model, const int iPart = -1);
  virtual void loop_by_particle(const Configuration& config, const ModelTwoBody& model, const int iPart = -1);
  double kloop_by_particle(const Configuration& config, const ModelTwoBody& model, const int iPart = -1) const;
  double kloop_by_particle(const Configuration& config, const ModelOneBody& model, const int iPart = -1) const;

  virtual void compute(const Configuration& config,
                      const ModelTwoBody& model,
                      const Select& selection);
  void compute(const Configuration& config,
              const ModelOneBody& model,
              const Select& selection);

  void init_loop() { energy_ = 0.; };
  /// evaluate pair interactions with a model

  int optimization_ = 1;

  double energy() const { return energy_; }
  void set_energy(const double energy) { energy_ = energy; }

  virtual ~VisitModel() {}

 private:
  double energy_;
  void benchmark_(const Configuration& config,
                  const ModelTwoBody& model,
                  const Select& selection);
};

}  // namespace feasst

#endif  // FEASST_CORE_VISIT_MODEL_H_
