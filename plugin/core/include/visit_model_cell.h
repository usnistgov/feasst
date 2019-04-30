
#ifndef FEASST_CORE_VISIT_MODEL_CELL_H_
#define FEASST_CORE_VISIT_MODEL_CELL_H_

#include "core/include/visit_model.h"

namespace feasst {

/**
  Compute many-body inter-particle interactions using a cell list.
  The cell index is determined by the number of times Configuration::init_cells
  is called.
 */
class VisitModelCell : public VisitModel {
 public:
  VisitModelCell() {}
  void compute(
      const ModelTwoBody& model,
      const ModelParams& model_params,
      Configuration * config,
      const int cell_index = 0) override;
  void compute(
      const ModelTwoBody& model,
      const ModelParams& model_params,
      const Select& selection,
      Configuration * config,
      const int cell_index = 0) override;
  std::shared_ptr<VisitModel> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  VisitModelCell(std::istream& istr);
  ~VisitModelCell() {}
 private:
  const std::string class_name_ = "VisitModelCell";
};

}  // namespace feasst

#endif  // FEASST_CORE_VISIT_MODEL_CELL_H_
