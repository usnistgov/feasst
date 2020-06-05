
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

class VisitModelInner {
 public:
  VisitModelInner() {}

  explicit VisitModelInner(const std::shared_ptr<EnergyMap> map) {
    set_energy_map(map); }

  virtual void compute(
    const int part1_index,
    const int site1_index,
    const int part2_index,
    const int site2_index,
    const Configuration * config,
    const ModelParams& model_params,
    const ModelTwoBody& model,
    const bool is_old_config,
    Position * relative,
    Position * pbc);

  virtual void precompute(Configuration * config) {
    if (energy_map_) {
      energy_map_->precompute(config);
    }
  }

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
      const Position * pbc) {
    energy_ += energy;
    if (energy_map_) {
      energy_map_->update(energy, part1_index, site1_index, site1_type,
        part2_index, site2_index, site2_type, squared_distance, pbc);
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
  void query_ixn(
      const int part1_index,
      const int site1_index,
      const int part2_index,
      const int site2_index) {
    if (energy_map_) {
      energy_ += energy_map_->query(part1_index, site1_index,
                                    part2_index, site2_index);
    }
  }

  double energy() const { return energy_; }

  void prep_for_revert(const Select& select) {
    if (energy_map_) {
      energy_map_->prep_for_revert(select);
    }
  }
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

  const EnergyMap& energy_map() const {
    return const_cast<EnergyMap&>(*energy_map_); }

  bool is_energy_map() const {
    if (energy_map_) { return true; } else { return false; } }

  bool is_energy_map_queryable() const;

  void check() const {
    if (energy_map_) {
      energy_map_->check();
    }
  }

  // serialize
  virtual void serialize(std::ostream& ostr) const {
    serialize_visit_model_inner_(ostr); }

  virtual std::shared_ptr<VisitModelInner> create(std::istream& istr) const {
    return std::make_shared<VisitModelInner>(istr); }
  std::map<std::string, std::shared_ptr<VisitModelInner> >& deserialize_map();
  std::shared_ptr<VisitModelInner> deserialize(std::istream& istr);
  explicit VisitModelInner(std::istream& istr);
  virtual ~VisitModelInner() {}

 protected:
  void serialize_visit_model_inner_(std::ostream& ostr) const;

 private:
  const std::string class_name_ = "VisitModelInner";
  double energy_ = 0.;
  double squared_distance_;
  std::shared_ptr<EnergyMap> energy_map_;
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
