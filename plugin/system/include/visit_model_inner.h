
#ifndef FEASST_SYSTEM_VISIT_MODEL_INNER_H_
#define FEASST_SYSTEM_VISIT_MODEL_INNER_H_

#include <map>
#include <string>
#include <memory>

namespace feasst {

class Configuration;
class EnergyMap;
class ModelParams;
class ModelTwoBody;
class ModelThreeBody;
class Position;
class Select;

typedef std::map<std::string, std::string> argtype;

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

  virtual void compute3body(
    const int part1_index,
    const int site1_index,
    const int part2_index,
    const int site2_index,
    const int part3_index,
    const int site3_index,
    const Position& r12,
    const Position& r13,
    const Configuration * config,
    const ModelParams& model_params,
    ModelThreeBody * model,
    const bool is_old_config,
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
    const Configuration& config);

  void clear_ixn(
      const int part1_index,
      const int site1_index,
      const int part2_index,
      const int site2_index);

  void query_ixn(const Select& select);

  double energy() const { return energy_; }

  void revert(const Select& select);

  void finalize(const Select& select);

  void set_energy_map(std::shared_ptr<EnergyMap> map);

  const EnergyMap& energy_map() const;

  bool is_energy_map() const;

  bool is_energy_map_queryable() const;

  void check(const Configuration& config) const;

  void synchronize_(const VisitModelInner& inner, const Select& perturbed);

  int cutoff_index() const { return cutoff_index_; }

  void set_skip_particle(const bool skip) { skip_particle_ = skip; }
  bool skip_particle() const { return skip_particle_; }

  /// Return true if update_interaction was called after compute
  int interacted() const { return interacted_; }
  void set_interacted(const int ixn) { interacted_ = ixn; }

  // serialize
  std::string class_name() const { return class_name_; }
  virtual void serialize(std::ostream& ostr) const;
  virtual std::shared_ptr<VisitModelInner> create(std::istream& istr) const {
    return std::make_shared<VisitModelInner>(istr); }
  virtual std::shared_ptr<VisitModelInner> create(argtype * args) const {
    return std::make_shared<VisitModelInner>(args); }
  std::map<std::string, std::shared_ptr<VisitModelInner> >& deserialize_map();
  std::shared_ptr<VisitModelInner> deserialize(std::istream& istr);
  std::shared_ptr<VisitModelInner> factory(const std::string name,
                                           argtype * args);
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
  int interacted_;

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
