
#ifndef FEASST_SYSTEM_VISIT_MODEL_INNER_H_
#define FEASST_SYSTEM_VISIT_MODEL_INNER_H_

#include <map>
#include <string>
#include <memory>
#include "system/include/model.h"
#include "system/include/energy_map.h"

namespace feasst {

class ModelParams;
class Configuration;
class ModelTwoBody;

// HWH rename to VisitInner (consider it can be used by BondVisitor, etc)
class VisitModelInner {
 public:
  explicit VisitModelInner(argtype args = argtype());
  explicit VisitModelInner(argtype * args);

  explicit VisitModelInner(const std::shared_ptr<EnergyMap> map) {
    set_energy_map(map); }

  virtual void compute(
    const int part1_index,
    const int site1_index,
    const int part2_index,
    const int site2_index,
    const Configuration * config,
    const ModelParams& model_params,
    ModelTwoBody * model,
    const bool is_old_config,
    Position * relative,
    Position * pbc,
    const double weight = 1.);

  virtual void precompute(Configuration * config);
  void set_energy(const double energy) { energy_ = energy; }
  void update_ixn(
      const double energy,
      const int part1_index,
      const int site1_index,
      const int site1_type,
      const int part2_index,
      const int site2_index,
      const int site2_type,
      const double squared_distance,
      const Position * pbc,
      const bool is_old_config,
      const Configuration& config) {
    energy_ += energy;
    if (energy_map_ && !is_old_config) {
      energy_map_->update(energy, part1_index, site1_index, site1_type,
        part2_index, site2_index, site2_type, squared_distance, pbc, config);
    }
  }
  void clear_ixn(
      const int part1_index,
      const int site1_index,
      const int part2_index,
      const int site2_index) {
    if (energy_map_) {
      energy_map_->clear(part1_index, site1_index, part2_index, site2_index);
    }
  }
  void query_ixn(const Select& select);

  double energy() const { return energy_; }

  void revert(const Select& select) {
    // HWH optimize, maybe map_new doens't have to be same
    // or have to revert, but how to calc new clusters
    // before finalize to check cluster constraint?
    if (energy_map_) {
      energy_map_->revert(select);
    }
  }
  void finalize(const Select& select) {
    if (energy_map_) {
      energy_map_->finalize(select);
    }
  }

  void set_energy_map(std::shared_ptr<EnergyMap> map) { energy_map_ = map; }

  const EnergyMap& energy_map() const;

  bool is_energy_map() const {
    if (energy_map_) { return true; } else { return false; } }

  bool is_energy_map_queryable() const;

  void check(const Configuration& config) const {
    if (energy_map_) {
      energy_map_->check(config);
    }
  }

  void synchronize_(const VisitModelInner& inner, const Select& perturbed) {
    if (energy_map_) {
      energy_map_->synchronize_(inner.energy_map(), perturbed);
    }
  }

  int cutoff_index() const { return cutoff_index_; }

  void set_skip_particle(const bool skip) { skip_particle_ = skip; }
  bool skip_particle() const { return skip_particle_; }

  // serialize
  std::string class_name() const { return class_name_; }
  virtual void serialize(std::ostream& ostr) const;
  virtual std::shared_ptr<VisitModelInner> create(std::istream& istr) const {
    return std::make_shared<VisitModelInner>(istr); }
  virtual std::shared_ptr<VisitModelInner> create(argtype * args) const {
    return std::make_shared<VisitModelInner>(args); }
  std::map<std::string, std::shared_ptr<VisitModelInner> >& deserialize_map();
  std::shared_ptr<VisitModelInner> deserialize(std::istream& istr);
  std::shared_ptr<VisitModelInner> factory(const std::string name, argtype * args);
  explicit VisitModelInner(std::istream& istr);
  virtual ~VisitModelInner() {}

 protected:
  std::string class_name_ = "VisitModelInner";
  void serialize_visit_model_inner_(std::ostream& ostr) const;

 private:
  double energy_ = 0.;
  double squared_distance_;
  int cutoff_index_ = -1;
  int cutoff_outer_index_ = -1;
  std::shared_ptr<EnergyMap> energy_map_;

  // temporariy and not serialized
  bool skip_particle_ = false;
};

inline std::shared_ptr<VisitModelInner> MakeVisitModelInner() {
  return std::make_shared<VisitModelInner>();
}

inline std::shared_ptr<VisitModelInner> MakeVisitModelInner(
    std::shared_ptr<EnergyMap> map) {
  return std::make_shared<VisitModelInner>(map);
}

}  // namespace feasst

#endif  // FEASST_SYSTEM_VISIT_MODEL_INNER_H_
