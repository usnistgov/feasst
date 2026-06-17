
#ifndef FEASST_CONFINEMENT_MODEL_TABLE_CARTESIAN_2D_INTEGR_H_
#define FEASST_CONFINEMENT_MODEL_TABLE_CARTESIAN_2D_INTEGR_H_

#include <vector>
#include "math/include/table.h"
#include "system/include/model_one_body.h"

namespace feasst {

class Domain;
class Random;
class Shape;
class System;

typedef std::map<std::string, std::string> argtype;

/**
  A tabular potential based on cartesian coordinates.
  Assumes symmetry along the x, y planes and that the Domain has no tilt.
  Integration of material does not take periodicity into account.
  E.g., the shapes extend forever and are not periodic in the domain.
 */
class ModelTableCart2DIntegr : public ModelOneBody {
 public:
  // Constructor for single site type tables.
  ModelTableCart2DIntegr(std::shared_ptr<Table2D> table);

  // Constructor for multiple site type tables.
  ModelTableCart2DIntegr(std::vector<std::shared_ptr<Table2D> > tables);

  const Table2D& table(const int site_type = 0) const;

  /// Generate the table by integration of the shape of the confinement over
  /// the entire and domain.
  void compute_table(
    Shape * shape,
    Domain * domain,
    Random * random,
    /// See Shape for documentation of integration_args.
    argtype integration_args,
    const int site_type = 0);

  /// Same as above, but parallelize the task with OMP.
  void compute_table_omp(
    Shape * shape,
    Domain * domain,
    Random * random,
    argtype integration_args,
    const int site_type = 0,
    /// See Thread for documentation of these two arguments.
    const int node = 0,
    const int num_node = 1);

  double energy(
    const Position& wrapped_site,
    const Site& site,
    const Configuration& config,
    const ModelParams& model_params) override;

  void serialize(std::ostream& ostr) const override;
  std::shared_ptr<Model> create(std::istream& istr) const override {
    return std::make_shared<ModelTableCart2DIntegr>(istr); }
  explicit ModelTableCart2DIntegr(std::istream& istr);
  virtual ~ModelTableCart2DIntegr() {}

 private:
  std::vector<std::shared_ptr<Table2D> > tables_;
};

inline std::shared_ptr<ModelTableCart2DIntegr> MakeModelTableCart2DIntegr(
    std::shared_ptr<Table2D> table) {
  return std::make_shared<ModelTableCart2DIntegr>(table);
}

}  // namespace feasst

#endif  // FEASST_CONFINEMENT_MODEL_TABLE_CARTESIAN_2D_INTEGR_H_
