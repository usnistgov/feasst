
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
  //@{
  /** @name Arguments
    - min_length: build cell list with given minimum distance between cells
      as a floating point number.
      If "max_sigma" then use the largest mixed sigma value.
      If "max_cutoff" then use the largest mixed cutoff value.
    - cell_group_index: compute cells only in given group index (default: 0).
    - cell_group: as above, but use the name of the group, not the index.
      Do not use at the same time as cell_group_index (default: "").
    - VisitModel arguments.
   */
  explicit VisitModelCell(argtype args);
  explicit VisitModelCell(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  /// Same as above, but with an inner.
  explicit VisitModelCell(std::shared_ptr<VisitModelInner> inner,
    argtype args);

  /// Return the cells.
  const Cells& cells() const { return cells_; }

  /// Return the unique cell number for the position.
  int cell_id(const Domain& domain, const Position& position) const;

  /// Same as above, but optimized.
  int cell_id_opt_(const Domain& domain, const Position& position);

  /// Same as base class, but also prepare the cells.
  void precompute(Configuration * config) override;

  void change_volume(const double delta_volume, const int dimension, Configuration * config) override;

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

  std::shared_ptr<VisitModel> create(std::istream& istr) const override {
    return std::make_shared<VisitModelCell>(istr); }
  std::shared_ptr<VisitModel> create(argtype * args) const override {
    return std::make_shared<VisitModelCell>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit VisitModelCell(std::istream& istr);
  VisitModelCell() {} // for mapper only
  virtual ~VisitModelCell() {}

  //@}
 private:
  Cells cells_;
  std::string min_length_;
  int group_index_;
  std::string group_;
  Position opt_origin_, opt_rel_, opt_pbc_;

  // temporary and not serialized
  Select one_site_select_;
  double opt_r2_;

  void position_tracker_(const Select& select, Configuration * config);
  double min_len_(const Configuration& config) const;
  void rebuild_(const Configuration& config);
};

inline std::shared_ptr<VisitModelCell> MakeVisitModelCell(
    argtype args = argtype()) {
  return std::make_shared<VisitModelCell>(args);
}

inline std::shared_ptr<VisitModelCell> MakeVisitModelCell(
    std::shared_ptr<VisitModelInner> inner,
    argtype args = argtype()) {
  return std::make_shared<VisitModelCell>(inner, args);
}

}  // namespace feasst

#endif  // FEASST_SYSTEM_VISIT_MODEL_CELL_H_
