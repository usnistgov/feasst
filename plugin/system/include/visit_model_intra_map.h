
#ifndef FEASST_SYSTEM_VISIT_MODEL_INTRA_MAP_H_
#define FEASST_SYSTEM_VISIT_MODEL_INTRA_MAP_H_

#include <map>
#include <string>
#include <vector>
#include <memory>
#include "system/include/visit_model.h"

namespace feasst {

typedef std::map<std::string, std::string> argtype;

/**
  Intra-particle interactions are computed here.
  As described in VisitModelIntra, this does not include bonded interaction
  energies, and only includes "inter"-like interactions.

  In this implementation, a map for each particle type is precomputed based
  on whether the "inter"-like interaction should be excluded due to the presence
  of a bond, angle or dihedral.
  In addition, dihedral interactions may be weighted.
  By default, all interactions are included except with self and there is no
  dihedral weight, so the user must specific if bonds, angles or dihedrals are
  to be excluded.
 */
class VisitModelIntraMap : public VisitModel {
 public:
  //@{
  /** @name Arguments
    - exclude_bonds: if true, exclude intra interactions between bonded sites
      (default: false).
    - exclude_angles: if true, exclude intra interactions between the two
      extremes of the angle sites (e.g, exclude AC of <ABC) (default: false).
    - exclude_dihedrals: if true, exclude intra interactions between the two
      extremes of the dihedral sites (e.g, exclude AD of <ABCD) (default: false).
    - dihedral_weight: if > 0, multiply intramolecular interactions for the
      1-4 dihedral sites by this weight (default: -1).
      For example, CHARMM and OPLS may set this weight to 1/2.
      If exclude_dihedral is true, this weight cannot be > 0.
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
  bool exclude_dihedrals_;
  double dihedral_weight_;
  std::vector<std::vector<std::vector<int> > > include_map_;
};

inline std::shared_ptr<VisitModelIntraMap> MakeVisitModelIntraMap(
    argtype args = argtype()) {
  return std::make_shared<VisitModelIntraMap>(args);
}

}  // namespace feasst

#endif  // FEASST_SYSTEM_VISIT_MODEL_INTRA_MAP_H_
