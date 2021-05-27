
#ifndef FEASST_SYSTEM_VISIT_MODEL_CELL_H_
#define FEASST_SYSTEM_VISIT_MODEL_CELL_H_

#include <memory>
#include "utils/include/arguments.h"
#include "system/include/visit_model.h"
#include "system/include/cells.h"

namespace feasst {

/**
  Compute many-body inter-particle interactions using a cell list.
 */
class VisitModelCell : public VisitModel {
 public:
  /**
    args:
    - min_length: build cell list with given minimum distance between cells.
    - cell_group: compute cells only in given group index (default: 0).
   */
  VisitModelCell(argtype args = argtype());

  /// Return the cells.
  const Cells& cells() const { return cells_; }

  /// Return the unique cell number for the position.
  int cell_id(const Domain& domain, const Position& position) const;

  /// Same as above, but optimized.
  int cell_id_opt_(const Domain& domain, const Position& position);

  void precompute(Configuration * config) override;

  void compute(
      ModelTwoBody * model,
      const ModelParams& model_params,
      Configuration * config,
      const int group_index) override;
  void compute(
      ModelTwoBody * model,
      const ModelParams& model_params,
      const Select& selection,
      Configuration * config,
      const int group_index) override;

  void finalize(const Select& select, Configuration * config) override;

  void check(const Configuration& config) const override;

  std::shared_ptr<VisitModel> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit VisitModelCell(std::istream& istr);
  virtual ~VisitModelCell() {}

 private:
  Cells cells_;
  double min_length_;
  int group_index_;
  Position opt_origin_, opt_rel_, opt_pbc_;

  // temporary and not serialized
  Select one_site_select_;
  double opt_r2_;

  void position_tracker_(const Select& select, Configuration * config);
};

inline std::shared_ptr<VisitModelCell> MakeVisitModelCell(
    argtype args = argtype()) {
  return std::make_shared<VisitModelCell>(args);
}

}  // namespace feasst

#endif  // FEASST_SYSTEM_VISIT_MODEL_CELL_H_
