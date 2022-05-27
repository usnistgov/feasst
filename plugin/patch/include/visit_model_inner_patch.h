
#ifndef FEASST_PATCH_VISIT_MODEL_INNER_PATCH_H_
#define FEASST_PATCH_VISIT_MODEL_INNER_PATCH_H_

#include "utils/include/arguments.h"
#include "system/include/visit_model.h"
#include "patch/include/patch_angle.h"

namespace feasst {

/**
  Patch interactions are defined by two sites. These are called centers or
  directors.

  The center site is used in the computation of the separation vector.
  Thus, for a Kern-Frenkel model, the centers would be the hard spheres.

  The director site is used to determine the orientation of the patch.

  The visitor should loop over a group_index which contains center sites only.
  Each pair of center sites is then subject to a cutoff which should be the
  maximum cutoff for all possible pairs of director types between the two
  centers.
  Once two centers are said to be within the cutoff, then all pair director
  interactions are computed.
  Directors are identified as those sites which are bonded to the center and
  whose bond is identified as a director bond.

  The orientational component of the patch is determined by this visitor class.
  But the distance-based component is determined by the model.
  The types for the model computation are determined by the directors.
  E.g., the sigmas, epsilons, cutoffs etc are set by the directors.

  The patch angle is defined by the angle between the center of mass separation
  vector of the non-director, center site, and the director.
  Thus, a 180 degree patch would encompass the entire surface.
 */
class VisitModelInnerPatch : public VisitModelInner {
 public:
  explicit VisitModelInnerPatch(argtype args = argtype());
  explicit VisitModelInnerPatch(argtype * args);
  void precompute(Configuration * config) override;
  void compute(
    const int part1_index,
    const int site1_index,
    const int part2_index,
    const int site2_index,
    const Configuration * config,
    const ModelParams& model_params,
    ModelTwoBody * model,
    const bool is_old_config,
    Position * relative,
    Position * pbc) override;

  const CosPatchAngle& cos_patch_angle() const { return cos_patch_angle_; }

  std::shared_ptr<VisitModelInner> create(std::istream& istr) const override {
    return std::make_shared<VisitModelInnerPatch>(istr); }
  std::shared_ptr<VisitModelInner> create(argtype * args) const override {
    return std::make_shared<VisitModelInnerPatch>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit VisitModelInnerPatch(std::istream& istr);
  virtual ~VisitModelInnerPatch() {}

 private:
  //int cos_patch_angle_index_;
  CosPatchAngle cos_patch_angle_;
  int director_index_ = -1;
};

inline std::shared_ptr<VisitModelInnerPatch> MakeVisitModelInnerPatch(
    argtype args = argtype()) {
  return std::make_shared<VisitModelInnerPatch>(args);
}

}  // namespace feasst

#endif  // FEASST_PATCH_VISIT_MODEL_INNER_PATCH_H_
