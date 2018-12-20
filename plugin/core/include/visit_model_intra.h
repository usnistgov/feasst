
#ifndef FEASST_CORE_VISIT_MODEL_INTRA_H_
#define FEASST_CORE_VISIT_MODEL_INTRA_H_

#include "core/include/visit_model.h"

namespace feasst {

/**
  Intra-particle interactions are computed here.
  this does not include bonded interaction energies, but "inter"-like models
  such as lennard jones but between sites in the same particle (e.g., long
  chains).
  Note that there is currently a simple method to determine which intraparticle
  interactions to consider.
  The intra_cut variable basically means to ignore any site index which is
  +/-intra_cut away.
  For example, by default, the intra_cut=1 is used to model a freely jointed
  linear chain
 */
class VisitModelIntra : public VisitModel {
 public:
  VisitModelIntra() { set_intra_cut(); }
  int intra_cut() const { return intra_cut_; }
  void set_intra_cut(const int cut = 1) { intra_cut_ = cut; }
  void compute(const Configuration& config,
      const ModelTwoBody& model,
      const Select& selection,
      const int group_index) override;
  ~VisitModelIntra() {}
 private:
  int intra_cut_;
};

}  // namespace feasst

#endif  // FEASST_CORE_VISIT_MODEL_INTRA_H_
