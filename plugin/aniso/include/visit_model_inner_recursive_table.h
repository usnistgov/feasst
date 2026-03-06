
#ifndef FEASST_ANISO_VISIT_MODEL_INNER_RECURSIVE_TABLE_H_
#define FEASST_ANISO_VISIT_MODEL_INNER_RECURSIVE_TABLE_H_

#include <string>
#include <sstream>
#include <memory>
#include "math/include/recursive_table.h"
#include "aniso/include/visit_model_inner_table.h"

namespace feasst {

class Configuration;

typedef std::map<std::string, std::string> argtype;

/**
 * This class is currently in development.
 */
class VisitModelInnerRecursiveTable : public VisitModelInnerTable {
 public:
  //@{
  /** @name Arguments
   */
  explicit VisitModelInnerRecursiveTable(argtype args = argtype());
  explicit VisitModelInnerRecursiveTable(argtype * args);

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
  double compute_aniso(const int type1, const int type2,
    const double squared_distance, const double s1, const double s2,
    const Configuration& config) const override;

  std::shared_ptr<VisitModelInner> create(std::istream& istr) const override {
    return std::make_shared<VisitModelInnerRecursiveTable>(istr); }
  std::shared_ptr<VisitModelInner> create(argtype * args) const override {
    return std::make_shared<VisitModelInnerRecursiveTable>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit VisitModelInnerRecursiveTable(std::istream& istr);
  virtual ~VisitModelInnerRecursiveTable() {}

  std::vector<std::vector<RecursiveTable5D> > contact_;
  std::vector<std::vector<RecursiveTable2D> > contact2d_;

  //@}
 private:
  std::string input_file_;
};

}  // namespace feasst

#endif  // FEASST_ANISO_VISIT_MODEL_INNER_RECURSIVE_TABLE_H_
