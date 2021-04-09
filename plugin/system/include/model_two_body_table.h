
#ifndef FEASST_SYSTEM_MODEL_TWO_BODY_TABLE_H_
#define FEASST_SYSTEM_MODEL_TWO_BODY_TABLE_H_

#include "utils/include/arguments.h"
#include "configuration/include/model_params.h"
#include "system/include/model_two_body.h"

namespace feasst {

class Table1D;

/**
  Tabulate two-body models and interpolate their interactions during the
  simulation.
 */
class ModelTwoBodyTable : public ModelTwoBody {
 public:
  /**
    args:
    - hard_sphere_threshold: when r < threshold, return NEAR_INFINITY
      (default: 0.2).
   */
  ModelTwoBodyTable(argtype args = argtype());
  ModelTwoBodyTable(argtype * args);

  /// Resize the table based on the number of site types
  void resize(const int num_site_types);

  /// Set a table potential between sites of type1 and type2.
  /// Table distances are normalized between 0 and 1, inclusive,
  /// where z = (r-rh)/(rc-rh), rh is the hard_sphere_threshold
  /// and rc is the cutoff.
  /// Thus, the first element of the table is for a distance of rh,
  /// and the last element of the table is for a distance of rc.
  void set(std::shared_ptr<Table1D> table,
    const int type1 = 0,
    const int type2 = 0);

  /// Tabulate an existing model.
  void set(const ModelParams& model_params,
    const int size,  /// size of table
    const int num_types,
    Model * model);

  /// Same as above, but with a shared_ptr model.
  void set(const ModelParams& model_params,
    const int size,  /// size of table
    const int num_types,
    std::shared_ptr<Model> model);

  /// Return the tabular potential.
  const Table1D& table(const int type1, const int type2) const {
    INFO(table_.size());
    INFO(table_[type1].size());
    return *table_[type1][type2]; }

  double energy(
    const double squared_distance,
    const int type1,
    const int type2,
    const ModelParams& model_params) override;

  std::shared_ptr<Model> create(std::istream& istr) const override {
    return std::make_shared<ModelTwoBodyTable>(istr); }
  void serialize(std::ostream& ostr) const override;
  explicit ModelTwoBodyTable(std::istream& istr);
  virtual ~ModelTwoBodyTable() {}

 protected:
  void serialize_model_two_body_table_(std::ostream& ostr) const;

 private:
  double hard_sphere_threshold_;
  std::vector<std::vector<std::shared_ptr<Table1D> > > table_;
};

inline std::shared_ptr<ModelTwoBodyTable> MakeModelTwoBodyTable(
  argtype args = argtype()) {
  return std::make_shared<ModelTwoBodyTable>(args);
}

}  // namespace feasst

#endif  // FEASST_SYSTEM_MODEL_TWO_BODY_TABLE_H_
