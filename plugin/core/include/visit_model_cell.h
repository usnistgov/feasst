
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
  void compute(
      const ModelTwoBody& model,
      Configuration * config,
      const int cell_index = 0) override;
  void compute(
      const ModelTwoBody& model,
      const Select& selection,
      Configuration * config,
      const int cell_index = 0) override;
  ~VisitModelCell() {}
 private:
};

}  // namespace feasst

#endif  // FEASST_CORE_VISIT_MODEL_CELL_H_
