
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
  VisitModelIntra() {}
  int intra_cut() const { return intra_cut_; }
  void set_intra_cut(const int cut) { intra_cut_ = cut; }
  void compute(
      const ModelTwoBody& model,
      const ModelParams& model_params,
      const Select& selection,
      Configuration * config,
      const int group_index) override;
  void compute(
      const ModelTwoBody& model,
      const ModelParams& model_params,
      Configuration * config,
      const int group_index) override {
    compute(model, model_params, config->selection_of_all(), config, group_index);
  }
  std::shared_ptr<VisitModel> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  VisitModelIntra(std::istream& istr);
  ~VisitModelIntra() {}

 private:
  const std::string class_name_ = "VisitModelIntra";
  int intra_cut_ = -1;
};

}  // namespace feasst

#endif  // FEASST_CORE_VISIT_MODEL_INTRA_H_
