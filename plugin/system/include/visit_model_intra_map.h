
#ifndef FEASST_SYSTEM_VISIT_MODEL_INTRA_MAP_H_
#define FEASST_SYSTEM_VISIT_MODEL_INTRA_MAP_H_

#include <vector>
#include <memory>
#include "utils/include/arguments.h"
#include "system/include/visit_model.h"

namespace feasst {

/**
  Intra-particle interactions are computed here.
  this does not include bonded interaction energies, but "inter"-like models
  such as lennard jones but between sites in the same particle (e.g., long
  chains).

  In this implementation, a map for each particle type is precomputed.
  For a given pair of site indices, the map returns true if the intra
  interaction is included.
 */
class VisitModelIntraMap : public VisitModel {
 public:
  //@{
  /** @name Arguments
   */

  /**
    By default, all interactions are included except with self.

    args:
    - exclude_bonds: if true, exclude intra interactions between bonded sites
      (default: false).
    - exclude_angles: if true, exclude intra interactions between the two
      extremes of the angle sites (e.g, exclude AC of <ABC) (default: false).
   */
  explicit VisitModelIntraMap(argtype args = argtype());
  explicit VisitModelIntraMap(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  /// Initialize include_map after VisitModel::precompute.
  void precompute(Configuration * config) override;

  /// Return 1 (true) if interactions between site1 and site2 in particle_type
  /// are included. Otherwise, return 0 (false);
  int include_map(const int particle_type, const int site1, const int site2) {
    return include_map_[particle_type][site1][site2]; }

  void compute(
      ModelTwoBody * model,
      const ModelParams& model_params,
      const Select& selection,
      Configuration * config,
      const int group_index) override;
  void compute(
      ModelTwoBody * model,
      const ModelParams& model_params,
      Configuration * config,
      const int group_index) override;
  std::shared_ptr<VisitModel> create(std::istream& istr) const override {
    return std::make_shared<VisitModelIntraMap>(istr); }
  std::shared_ptr<VisitModel> create(argtype * args) const override {
    return std::make_shared<VisitModelIntraMap>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit VisitModelIntraMap(std::istream& istr);
  ~VisitModelIntraMap() {}

  //@}
 private:
  bool exclude_bonds_;
  bool exclude_angles_;
  std::vector<std::vector<std::vector<int> > > include_map_;
};

inline std::shared_ptr<VisitModelIntraMap> MakeVisitModelIntraMap(
    argtype args = argtype()) {
  return std::make_shared<VisitModelIntraMap>(args);
}

}  // namespace feasst

#endif  // FEASST_SYSTEM_VISIT_MODEL_INTRA_MAP_H_
