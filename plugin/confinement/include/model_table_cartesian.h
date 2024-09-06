
#ifndef FEASST_CONFINEMENT_MODEL_TABLE_CARTESIAN_H_
#define FEASST_CONFINEMENT_MODEL_TABLE_CARTESIAN_H_

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
  A tabular potential based on cartesian coordinates.
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

/**
  A tabular potential based on cartesian coordinates.
  Assumes symmetry along the x, y and z planes and that the Domain has no tilt.
  Integration of material does not take periodicity into account.
  E.g., the shapes extend forever and are not periodic in the domain.
 */
class ModelTableCart3DIntegr : public ModelOneBody {
 public:
  // Constructor for single site type tables.
  explicit ModelTableCart3DIntegr(std::shared_ptr<Table3D> table);

  // Constructor for multiple site type tables.
  explicit ModelTableCart3DIntegr(std::vector<std::shared_ptr<Table3D> > tables);

  //@{
  /** @name Arguments
    - table_file: file name for the table.
    - shape_file: ShapeFile that describes the shape.
    - ModelTableCart3DIntegr::compute_table::integration_args.
    - use_omp: use OpenMP to compute the table (default: false).
    - node: for parallelization, see compute_table_omp (default: 0).
    - num_node: for parallelization, see compute_table_omp (default: 1).
    - Table3D arguments.
    - Domain arguments.

    The format for the table file is as follows.

    The first line should be 'site_types' followed by the number of site types
    and then the identity of each of those site types in order of the tables
    given below. (e.g., "site_types n i" where n is the number of site types and
    each following number is the type of each anisotropic site.)

    The remaining lines are the individual tables for each of the site types.
   */
  explicit ModelTableCart3DIntegr(argtype args = argtype());
  explicit ModelTableCart3DIntegr(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  /// Return the table for a given site type.
  const Table3D& table(const int site_type = 0) const;

  /**
    Generate the table by integration of a shape, which represents a continuous
    medium, over the entire domain.
   */
  void compute_table(
    Shape * shape,
    const Domain& domain,
    Random * random,
    /// See Shape for documentation of integration_args.
    argtype integration_args,
    const int site_type = 0);

  /// Same as above, but parallelize the task with OMP.
  void compute_table_omp(
    Shape * shape,
    const Domain& domain,
    Random * random,
    argtype integration_args,
    const int site_type = 0,
    /// See Thread for documentation of these two arguments.
    const int node = 0,
    const int num_nodes = 1);

  // Same as above, but while parsing arguments
  void compute_table(
    Shape * shape,
    const Domain& domain,
    Random * random,
    /// See Shape for documentation of integration_args.
    argtype * integration_args,
    const int site_type = 0);

  /// Same as above, but parallelize the task with OMP.
  void compute_table_omp(
    Shape * shape,
    const Domain& domain,
    Random * random,
    argtype * integration_args,
    const int site_type = 0,
    /// See Thread for documentation of these two arguments.
    const int node = 0,
    const int num_nodes = 1);

  // HWH refactor so that site_types are harvested from select
  /**
    Generate the table by computing the energy of interaction of the select
    with the rest of the system.
    The select is assumed to be a single site, so that tables can be generated
    for each site type.
   */
  void compute_table(
    System * system,
    Select * select,
    const int site_type = 0);

  /// Same as above, but parallelize the task with OMP
  void compute_table_omp(
    System * system,
    Select * select,
    const int site_type = 0,
    /// See Thread for documentation of these two arguments.
    const int node = 0,
    const int num_node = 1);

  double energy(
    const Position& wrapped_site,
    const Site& site,
    const Configuration& config,
    const ModelParams& model_params) override;

  void write(const std::string file_name) const;
  void read(const std::string file_name);

  void serialize(std::ostream& ostr) const override;
  std::shared_ptr<Model> create(std::istream& istr) const override {
    return std::make_shared<ModelTableCart3DIntegr>(istr); }
  explicit ModelTableCart3DIntegr(std::istream& istr);
  virtual ~ModelTableCart3DIntegr() {}

 private:
  std::vector<std::shared_ptr<Table3D> > tables_;
  std::vector<int> site_types_;
  argtype args_;
};

inline std::shared_ptr<ModelTableCart3DIntegr> MakeModelTableCart3DIntegr(
    std::shared_ptr<Table3D> table) {
  return std::make_shared<ModelTableCart3DIntegr>(table);
}

inline std::shared_ptr<ModelTableCart3DIntegr> MakeModelTableCart3DIntegr(
    argtype args = argtype()) {
  return std::make_shared<ModelTableCart3DIntegr>(args);
}

}  // namespace feasst

#endif  // FEASST_CONFINEMENT_MODEL_TABLE_CARTESIAN_H_
