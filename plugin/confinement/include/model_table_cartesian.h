
#ifndef FEASST_CONFINEMENT_MODEL_TABLE_CARTESIAN_H_
#define FEASST_CONFINEMENT_MODEL_TABLE_CARTESIAN_H_

#include <vector>
#include "utils/include/arguments.h"
#include "math/include/table.h"
#include "system/include/model_one_body.h"

namespace feasst {

// HWH have domain return scaled coordinates?
/**
  A tabular potential based on cartesian coordinates.
  Assumes symmetry along the x, y and z planes and that the Domain has no tilt.
 */
class ModelTableCart3FoldSym : public ModelOneBody {
 public:
  ModelTableCart3FoldSym(std::shared_ptr<Table3D> table) { table_ = table; }
  double energy(
      const Site& site,
      const Configuration& config,
      const ModelParams& model_params) const override;
  void serialize(std::ostream& ostr) const override;
  std::shared_ptr<Model> create(std::istream& istr) const override {
    return std::make_shared<ModelTableCart3FoldSym>(istr); }
  explicit ModelTableCart3FoldSym(std::istream& istr);
  virtual ~ModelTableCart3FoldSym() {}

 private:
  const std::string class_name_ = "ModelTableCart3FoldSym";
  std::shared_ptr<Table3D> table_;
};

inline std::shared_ptr<ModelTableCart3FoldSym> MakeModelTableCart3FoldSym(
    std::shared_ptr<Table3D> table) {
  return std::make_shared<ModelTableCart3FoldSym>(table);
}

}  // namespace feasst

#endif  // FEASST_CONFINEMENT_MODEL_TABLE_CARTESIAN_H_
