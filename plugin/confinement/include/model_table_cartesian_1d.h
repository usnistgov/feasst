
#ifndef FEASST_CONFINEMENT_MODEL_TABLE_CARTESIAN_1D_H_
#define FEASST_CONFINEMENT_MODEL_TABLE_CARTESIAN_1D_H_

#include <vector>
#include "math/include/table.h"
#include "system/include/model_one_body.h"

namespace feasst {

typedef std::map<std::string, std::string> argtype;

void precompute_table1D(const Configuration& config, std::string table_file,
  std::vector<std::unique_ptr<Table1D> > * tables);

// HWH: Add HalfSpace or some other Shape instead of dimension argument,
// or maybe make it work with all Shapes? but how to scale the coordinate?

/**
  A tabular potential for interactions in one Cartesian dimension.
  The table file is formatted as follows:

  The first line should be "site_types=" followed by a comma-separated list of
  the labels of each site type.
  The remaining lines are the individual tables for each of the site types.
  For example, "site_types=O,H" will then contain a line for the table of "O"
  sites, then another line for the table of "H" sites.
  Each table is given in a single line by a number of space-separated values
  that represent the interaction energy spanning the range of a coordinate, z,
  scaled by the periodic boundary lengths, L, with the range z=[-0.5, 0.5] and
  linear spacing.
  For example, if the line contains 5 values, then z=[-0.5, -0.25, 0, 0.25, 0.5]
  and the potential will be interpolated with
  Table1D::forward_difference_interpolation.
 */
class ModelTableCart1D : public ModelOneBody {
 public:
  //@{
  /** @name Arguments
    - table_file: file name for the table.
    - dimension: the cartesian dimension to apply the Potential (default: 2).
   */
  explicit ModelTableCart1D(argtype args = argtype());
  explicit ModelTableCart1D(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  void precompute(Configuration * config, ModelParams * params) override {
    precompute_table1D(*config, table_file_, &tables_); }

  double energy(
    const Position& wrapped_site,
    const Site& site,
    const Configuration& config,
    const ModelParams& model_params) override;

  void serialize(std::ostream& ostr) const override;
  std::shared_ptr<Model> create(std::istream& istr) const override {
    return std::make_shared<ModelTableCart1D>(istr); }
  std::shared_ptr<Model> create(argtype * args) const override {
    return std::make_shared<ModelTableCart1D>(args); }
  explicit ModelTableCart1D(std::istream& istr);
  virtual ~ModelTableCart1D();

  //@}
 private:
  std::string table_file_;
  int dimension_;
  std::vector<std::unique_ptr<Table1D> > tables_;
};

}  // namespace feasst

#endif  // FEASST_CONFINEMENT_MODEL_TABLE_CARTESIAN_1D_H_
