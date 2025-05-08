
#ifndef FEASST_SYSTEM_VISIT_MODEL_H_
#define FEASST_SYSTEM_VISIT_MODEL_H_

#include <memory>
#include <string>
#include <map>
#include "math/include/position.h"
#include "system/include/synchronize_data.h"

namespace feasst {

class Domain;
class Configuration;
class Model;
class ModelParams;
class ModelOneBody;
class ModelTwoBody;
class ModelThreeBody;
class Select;
class VisitModelInner;

struct pair3body {
  int site1, site2, part1, part2;
  Position rel;
};

typedef std::map<std::string, std::string> argtype;

/**
  This class loops a Model over a Configuration or Select.

  See Model for a description of the compute methods. These are mirrored by
  simply switching the calling object and the first argument
  (.e.g, Model.compute(Visitor, ...) vs Visitor.compute(Model, ...)
 */
class VisitModel {
 public:
  VisitModel();  // use the default VisitModelInner
  explicit VisitModel(std::shared_ptr<VisitModelInner> inner);

  //@{
  /** @name Arguments
    - energy_cutoff: energy above this value will immediately end loop without
      computing the energy of the remaining sites in the loop.
      Must be > 1e10 because too low could result in an accepted trial.
      If -1, ignore energy_cutoff (default: -1).
    - VisitModelInner: derived class VisitModelInner (default: VisitModelInner).
   */
  explicit VisitModel(argtype args);
  explicit VisitModel(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  double energy_cutoff() const { return energy_cutoff_; }

  void set_inner(const std::shared_ptr<VisitModelInner> inner);

  const VisitModelInner& inner() const;

  virtual void compute(
      ModelOneBody * model,
      const ModelParams& model_params,
      Configuration * config,
      const int group_index = 0);
  virtual void compute(
      ModelOneBody * model,
      const ModelParams& model_params,
      const Select& selection,
      Configuration * config,
      const int group_index = 0);
  virtual void compute(
      ModelTwoBody * model,
      const ModelParams& model_params,
      Configuration * config,
      const int group_index = 0);
  virtual void compute(
      ModelTwoBody * model,
      const ModelParams& model_params,
      const Select& selection,
      Configuration * config,
      const int group_index = 0);
  virtual void compute(
      ModelThreeBody * model,
      const ModelParams& model_params,
      Configuration * config,
      const int group_index = 0);
  virtual void compute(
      ModelThreeBody * model,
      const ModelParams& model_params,
      const Select& selection,
      Configuration * config,
      const int group_index = 0);

  // If model parameters are not given, then obtain them from the configuration.
  void compute(
      ModelOneBody * model,
      Configuration * config,
      const int group_index = 0);
  void compute(
      ModelOneBody * model,
      const Select& selection,
      Configuration * config,
      const int group_index = 0);
  virtual void compute(
      ModelTwoBody * model,
      Configuration * config,
      const int group_index = 0);
  void compute(
      ModelTwoBody * model,
      const Select& selection,
      Configuration * config,
      const int group_index = 0);
  void compute(
      ModelThreeBody * model,
      Configuration * config,
      const int group_index = 0);
  void compute(
      ModelThreeBody * model,
      const Select& selection,
      Configuration * config,
      const int group_index = 0);

  // compute interactions between particles in the selection
  void compute_between_selection(
    ModelTwoBody * model,
    const ModelParams& model_params,
    const Select& selection,
    Configuration * config,
    const bool is_old_config,
    Position * relative,
    Position * pbc);

  /// Return the energy.
  double energy() const { return energy_; }

  /// Set the energy.
  void set_energy(const double energy) { energy_ = energy; }

  void zero_energy();

  /// Increment the energy.
  void increment_energy(const double energy) { energy_ += energy; }

  /// Test if energy of whole system is consistent with sum of energy
  /// of selection by particles.
  void check_energy(
      Model * model,
      Configuration * config,
      const int group_index = 0);

  virtual void revert(const Select& select);
  virtual void finalize(const Select& select, Configuration * config);

  virtual void precompute(Configuration * config);

  virtual void check(const Configuration& config) const;

  // Synchronize with another object of the same type.
  // Typically used with prefetch.
  virtual void synchronize_(const VisitModel& visit, const Select& perturbed);
  const SynchronizeData& data() const { return data_; }
  const SynchronizeData& manual_data() const { return manual_data_; }

  /// Change the volume.
  virtual void change_volume(const double delta_volume, const int dimension,
    Configuration * config) {}

  /// Return the ModelParams index of epsilon.
  int epsilon_index() const { return epsilon_index_; }

  /// Return the ModelParams index of sigma.
  int sigma_index() const { return sigma_index_; }

  /// Return the ModelParams index of cutoff.
  int cutoff_index() const { return cutoff_index_; }

  /// Return the ModelParams index of charge.
  int charge_index() const { return charge_index_; }

  // serialization
  std::string class_name() const { return class_name_; }
  virtual void serialize(std::ostream& ostr) const;
  virtual std::shared_ptr<VisitModel> create(std::istream& istr) const {
    return std::make_shared<VisitModel>(istr); }
  virtual std::shared_ptr<VisitModel> create(argtype * args) const {
    return std::make_shared<VisitModel>(args); }
  std::map<std::string, std::shared_ptr<VisitModel> >& deserialize_map();
  std::shared_ptr<VisitModel> deserialize(std::istream& istr);
  std::shared_ptr<VisitModel> factory(const std::string name, argtype * args);
  explicit VisitModel(std::istream& istr);
  virtual ~VisitModel();

  //@}
 protected:
  std::string class_name_ = "VisitModel";
  void serialize_visit_model_(std::ostream& ostr) const;
  VisitModelInner * get_inner_() const;

  // HWH hacky addition for optimization: also, prep inner for reverting,
  // because this is called at beginning of every pair-wise selection compute
  // optimization to avoid repeated construction of Position.
  std::shared_ptr<Position> relative_, pbc_, origin_;
  void init_relative_(const Domain& domain);

  SynchronizeData data_;  // all data is copied at synchronization
  SynchronizeData manual_data_;  // data is manually copied

  // Three body list of pairs
  std::vector<pair3body> pairs_;
  void record_pair_(const int part1_index, const int site1_index, const int part2_index, const int site2_index, const Position& rel, int * num_pair, VisitModelInner * inner);
  void pair_pair_(const int num_pair,
    ModelThreeBody * model,
    const ModelParams& model_params,
    Configuration * config,
    const Select * sel);
  bool find_in_pair3body(const int ipart, const int isite, const int num_pair);
  bool is_old_config_(const Select& selection) const;

  // If possible, query energy map of old configuration instead of pair loop
  bool is_queryable_(const Select& selection, const bool is_old_config, VisitModelInner * inner);

 private:
  double energy_ = 0.;
  std::shared_ptr<VisitModelInner> inner_;
  int epsilon_index_ = -1;
  int sigma_index_ = -1;
  int cutoff_index_ = -1;
  int charge_index_ = -1;
  double energy_cutoff_;
};

inline std::shared_ptr<VisitModel> MakeVisitModel(argtype args = argtype()) {
  return std::make_shared<VisitModel>(args);
}

inline std::shared_ptr<VisitModel> MakeVisitModel(
    std::shared_ptr<VisitModelInner> inner) {
  return std::make_shared<VisitModel>(inner);
}

}  // namespace feasst

#endif  // FEASST_SYSTEM_VISIT_MODEL_H_
