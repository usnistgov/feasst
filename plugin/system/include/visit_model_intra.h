
#ifndef FEASST_SYSTEM_VISIT_MODEL_INTRA_H_
#define FEASST_SYSTEM_VISIT_MODEL_INTRA_H_

#include <memory>
#include "utils/include/arguments.h"
#include "system/include/visit_model.h"

namespace feasst {

/**
  Intra-particle interactions are computed here.
  this does not include bonded interaction energies, but "inter"-like models
  such as lennard jones but between sites in the same particle (e.g., long
  chains).
 */
class VisitModelIntra : public VisitModel {
 public:
  /**
    args:
    - intra_cut: ignore the interaction between a pair of sites when the difference
      between their indices, |i-j| <= intra_cut (integer, default: -1).
   */
  explicit VisitModelIntra(argtype args = argtype());
  explicit VisitModelIntra(argtype * args);
  int intra_cut() const { return intra_cut_; }
  void set_intra_cut(const int cut) { intra_cut_ = cut; }
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
    return std::make_shared<VisitModelIntra>(istr); }
  std::shared_ptr<VisitModel> create(argtype * args) const override {
    return std::make_shared<VisitModelIntra>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit VisitModelIntra(std::istream& istr);
  ~VisitModelIntra() {}

 private:
  int intra_cut_;
};

inline std::shared_ptr<VisitModelIntra> MakeVisitModelIntra(
    argtype args = argtype()) {
  return std::make_shared<VisitModelIntra>(args);
}

}  // namespace feasst

#endif  // FEASST_SYSTEM_VISIT_MODEL_INTRA_H_
