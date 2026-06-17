
#ifndef FEASST_CONFINEMENT_MODEL_TABLE_CARTESIAN_1D_HARD_H_
#define FEASST_CONFINEMENT_MODEL_TABLE_CARTESIAN_1D_HARD_H_

#include <vector>
#include "math/include/table.h"
#include "system/include/model_one_body.h"

namespace feasst {

typedef std::map<std::string, std::string> argtype;

class Domain;
class Random;
class Shape;
class System;

// HWH have domain return scaled coordinates?
/**
  A tabular potential for a hard surface based on x,y cartesian coordinates.
  Assumes symmetry along the x plane and that the Domain has no tilt.
 */
class ModelTableCart1DHard : public ModelOneBody {
 public:
  // Constructor for single site type tables.
  ModelTableCart1DHard(std::shared_ptr<Table1D> table);

  // Constructor for multiple site type tables.
  ModelTableCart1DHard(std::vector<std::shared_ptr<Table1D> > tables);

  const Table1D& table(const int site_type = 0) const;

  /**
    Generate the table by finding where the point is inside the shape and the
    nearest distance to the surface is half of the diameter.
    The initial bounds are [0, L/2] inclusive, assuming
    a plane (or line) of symmetry at origin perpendicular to y axis.

    args:
    - diameter: diameter of the sphere (default: 1)
   */
  void compute_table(
    Shape * shape,
    Domain * domain,
    Random * random,
    argtype args = argtype(),
    const int site_type = 0);

  double energy(
    const Position& wrapped_site,
    const Site& site,
    const Configuration& config,
    const ModelParams& model_params) override;

  void serialize(std::ostream& ostr) const override;
  std::shared_ptr<Model> create(std::istream& istr) const override {
    return std::make_shared<ModelTableCart1DHard>(istr); }
  explicit ModelTableCart1DHard(std::istream& istr);
  virtual ~ModelTableCart1DHard() {}

 private:
  std::vector<std::shared_ptr<Table1D> > tables_;
};

inline std::shared_ptr<ModelTableCart1DHard> MakeModelTableCart1DHard(
    std::shared_ptr<Table1D> table) {
  return std::make_shared<ModelTableCart1DHard>(table);
}

}  // namespace feasst

#endif  // FEASST_CONFINEMENT_MODEL_TABLE_CARTESIAN_1D_HARD_H_
