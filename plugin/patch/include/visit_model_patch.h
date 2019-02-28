
#ifndef FEASST_CORE_VISIT_MODEL_PATCH_H_
#define FEASST_CORE_VISIT_MODEL_PATCH_H_

#include "core/include/visit_model.h"

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
 */
class VisitModelPatch : public VisitModel {
 public:
  // HWH make a helper function which
  // 1. sets up rcut of centers based on patches
  // assigns group for centers and uses that
  // and finally, stores model param index for patch angles

  // add all centers that are bonded to directors as a group for group index
//  int centers_group(System * system) {
//    Configuration * config = system()->get_configuration();
//    const int group_index = config->num_groups();
//    Group group;
//    std::vector<int> types;
//    const ModelParams& params = config->model_params();
//    for (int site_type = 0; site_type < config->num_site_types(); ++site_type) {
//      params.
//      
//      if (
//      
//      config->add(Group()

//    }
//  }

  // HWH note: terrible interface to patch angle for all types
  double cpa_sq_ = pow(cos(90./180.*PI), 2);

 protected:
  // compute the interaction between a pair of centers
  void inner_(
      const int part1_index,
      const int site1_index,
      const int part2_index,
      const int site2_index,
      const Configuration * config,
      const ModelParams& model_params,
      const ModelTwoBody& model,
      Position * relative) override;
};

}  // namespace feasst

#endif  // FEASST_CORE_VISIT_MODEL_PATCH_H_
