
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

  /// Add selection to cells.
  void add_to_cell_list(const Select& select,
                        const int particle_cell) {
    cells_.add(select, particle_cell); }

  /// Update selection in cells.
  void update_cell_list(const Select& select,
                        const int cell_new,
                        const int cell_old) {
    cells_.update(select, cell_new, cell_old); }

  /// Remove selection from cells.
  void remove_from_cell_list(const Select& select,
                             const int cell) {
    cells_.remove(select, cell); }

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

  void finalize(const Select& select) override;

  void check() const override;

  std::shared_ptr<VisitModel> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit VisitModelCell(std::istream& istr);
  virtual ~VisitModelCell() {}

 private:
  Cells cells_;
  double min_length_;
  int group_index_;
  int cell_index_ = 0;  // HWH update for more than one cell list
  Position opt_origin_, opt_rel_, opt_pbc_;

  // temporary and not serialized
  Configuration * config_;
  Select one_site_select_;
  double opt_r2_;

  void position_tracker_(const Select& select);
};

inline std::shared_ptr<VisitModelCell> MakeVisitModelCell(
    argtype args = argtype()) {
  return std::make_shared<VisitModelCell>(args);
}

}  // namespace feasst

#endif  // FEASST_SYSTEM_VISIT_MODEL_CELL_H_
