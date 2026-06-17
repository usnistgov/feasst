
#ifndef FEASST_ANISO_MODEL_RECURSIVE_TABLE_H_
#define FEASST_ANISO_MODEL_RECURSIVE_TABLE_H_

#include <memory>
#include <vector>
#include "math/include/recursive_table.h"
#include "system/include/model_two_body.h"

namespace feasst {

//class RecursiveTable1D;
class Table1D;

typedef std::map<std::string, std::string> argtype;

/**
  This class is experimental.
  This can have tables inside of tables (RecursiveTable).
  It serves more as a demonstration than a recommended method, because
  alternative methods such as splines would probably be better.
 */
class ModelRecursiveTable : public ModelTwoBody {
 public:
  //@{
  /** @name Arguments
    - input_file: checkpoint file from BuildRecursiveTable.
      For more than one site type, provide a comma-separated list of the colon-
      separated pairs of "[site type name]_[site type name]:filename"
      For example: "A_A:file1,A_B:file2,B_B:file2" for A and B interactions.
   */
  explicit ModelRecursiveTable(argtype args = argtype());
  explicit ModelRecursiveTable(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  void precompute(Configuration * config) override;

  double energy(
    const double squared_distance,
    const int type1,
    const int type2,
    const ModelParams& model_params) override;

  std::shared_ptr<Model> create(std::istream& istr) const override {
    return std::make_shared<ModelRecursiveTable>(istr); }
  std::shared_ptr<Model> create(argtype * args) const override {
    return std::make_shared<ModelRecursiveTable>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit ModelRecursiveTable(std::istream& istr);
  virtual ~ModelRecursiveTable();

  //@}
  std::vector<std::vector<double> > lower_;
  std::vector<std::vector<double> > upper_;
  std::vector<int> t2index_;
  std::vector<std::vector<RecursiveTable1D> > energy_;

 private:
  std::string input_file_;
};

void files_to_types_recursive_(const std::vector<std::string>& files, const Configuration& config, std::vector<int> * t2index, std::vector<std::string> * site_type_names, std::vector<int> * site_types);
std::string filename_to_idx_recursive_(const std::string& filename, const Configuration& config, const std::vector<int>& t2index, int * idx1, int * idx2);

}  // namespace feasst

#endif  // FEASST_ANISO_MODEL_RECURSIVE_TABLE_H_
