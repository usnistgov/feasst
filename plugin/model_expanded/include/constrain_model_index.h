
#ifndef FEASST_MODEL_EXPANDED_CONSTRAIN_MODEL_INDEX_H_
#define FEASST_MODEL_EXPANDED_CONSTRAIN_MODEL_INDEX_H_

#include <memory>
#include "monte_carlo/include/constraint.h"

namespace feasst {

typedef std::map<std::string, std::string> argtype;

/**
  Constrain Model::model_index.
 */
class ConstrainModelIndex : public Constraint {
 public:
  //@{
  /** @name Arguments
    - maximum: maximum model_index. If -1, no limit (default: -1).
    - minimum: minimum model_index (default: 0).
    - potential_index: index of potential for model (default: 0).
   */
  explicit ConstrainModelIndex(argtype args = argtype());
  explicit ConstrainModelIndex(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  /// Return the potential_index.
  int potential_index() const { return potential_index_; }

  bool is_allowed(const System& system,
    const Criteria& criteria,
    const Acceptance& acceptance) const override;

  /// Return the Model::model_index proposed by a Trial.
  int model_index(const System& system, const Acceptance& acceptance) const;

  std::shared_ptr<Constraint> create(std::istream& istr) const override {
    return std::make_shared<ConstrainModelIndex>(istr); }
  std::shared_ptr<Constraint> create(argtype * args) const override {
    return std::make_shared<ConstrainModelIndex>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit ConstrainModelIndex(std::istream& istr);
  virtual ~ConstrainModelIndex() {}

  //@}
 protected:
  void serialize_constrain_model_index_(std::ostream& ostr) const;
  int maximum_;
  int minimum_;
  int potential_index_;
};

inline std::shared_ptr<ConstrainModelIndex> MakeConstrainModelIndex(
    argtype args = argtype()) {
  return std::make_shared<ConstrainModelIndex>(args);
}

}  // namespace feasst

#endif  // FEASST_MODEL_EXPANDED_CONSTRAIN_MODEL_INDEX_H_
