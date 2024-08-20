
#ifndef FEASST_ANISO_VISIT_MODEL_INNER_NN_H_
#define FEASST_ANISO_VISIT_MODEL_INNER_NN_H_

#include <memory>
#include "math/include/table.h"
#include "math/include/matrix.h"
#include "math/include/euler.h"
#include "aniso/include/visit_model_inner_table.h"

namespace feasst {

typedef std::map<std::string, std::string> argtype;

/**
 * This class is currently in development.
 */
class VisitModelInnerNN : public VisitModelInnerTable {
 public:
  //@{
  /** @name Arguments
   */
  explicit VisitModelInnerNN(argtype args = argtype());
  explicit VisitModelInnerNN(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  void precompute_cutoffs(Configuration * config) override {};
  void read_table(const std::string table_file,
    const bool ignore_energy, Configuration * config) override;
  double compute_aniso(const int type1, const int type2,
    const double squared_distance, const double s1, const double s2,
    const double e1, const double e2, const double e3, const Configuration& config) const override;

  std::shared_ptr<VisitModelInner> create(std::istream& istr) const override {
    return std::make_shared<VisitModelInnerNN>(istr); }
  std::shared_ptr<VisitModelInner> create(argtype * args) const override {
    return std::make_shared<VisitModelInnerNN>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit VisitModelInnerNN(std::istream& istr);
  virtual ~VisitModelInnerNN() {}

  //@}
 private:
  // HWH declare NN objects here.
};

inline std::shared_ptr<VisitModelInnerNN> MakeVisitModelInnerNN(
    argtype args = argtype()) {
  return std::make_shared<VisitModelInnerNN>(args);
}

}  // namespace feasst

#endif  // FEASST_ANISO_VISIT_MODEL_INNER_NN_H_
