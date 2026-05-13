
#ifndef FEASST_ANISO_RECURSIVE_TABLE_H_
#define FEASST_ANISO_RECURSIVE_TABLE_H_

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
class RecursiveTable : public VisitModelInnerTable {
 public:
  //@{
  /** @name Arguments
    - VisitModelInnerTable arguments.
   */
  explicit RecursiveTable(argtype args = argtype());
  explicit RecursiveTable(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  void precompute_cutoffs(Configuration * config) override {};
  void read_table(const std::string table_file,
    const bool ignore_energy, Configuration * config) override;
  double compute_aniso(const int type1, const int type2,
    const double squared_distance, const double s1, const double s2,
    const double e1, const double e2, const double e3,
    const Configuration& config, const ModelParams& model_params,
    const int stype1, const int stype2) const override;
  double compute_aniso(const int type1, const int type2,
    const double squared_distance, const double s1, const double s2,
    const Configuration& config, const ModelParams& model_params,
    const int stype1, const int stype2) const override;
  double compute_aniso(const int type1, const int type2,
    const double squared_distance, const double s1,
    const Configuration& config, const ModelParams& model_params,
    const int stype1, const int stype2) const override;

  std::shared_ptr<VisitModelInner> create(std::istream& istr) const override {
    return std::make_shared<RecursiveTable>(istr); }
  std::shared_ptr<VisitModelInner> create(argtype * args) const override {
    return std::make_shared<RecursiveTable>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit RecursiveTable(std::istream& istr);
  virtual ~RecursiveTable() {}

  std::vector<std::vector<RecursiveTable5D> > contact_;
  std::vector<std::vector<RecursiveTable5D> > cutoff_;
  std::vector<std::vector<RecursiveTable6D> > energy_;
  std::vector<std::vector<RecursiveTable2D> > contact2d_;
  std::vector<std::vector<RecursiveTable2D> > cutoff2d_;
  std::vector<std::vector<RecursiveTable3D> > energy3d_;
  std::vector<std::vector<RecursiveTable1D> > contact1d_;
  std::vector<std::vector<RecursiveTable1D> > cutoff1d_;
  std::vector<std::vector<RecursiveTable2D> > energy2d_;

  //@}
 private:
  std::string input_file_;
};

}  // namespace feasst

#endif  // FEASST_ANISO_RECURSIVE_TABLE_H_
