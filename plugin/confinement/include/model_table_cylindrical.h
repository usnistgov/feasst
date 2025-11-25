
#ifndef FEASST_CONFINEMENT_MODEL_TABLE_CYLINDRICAL_H_
#define FEASST_CONFINEMENT_MODEL_TABLE_CYLINDRICAL_H_

#include <memory>
#include <string>
#include <vector>
#include "system/include/model_one_body.h"

namespace feasst {

typedef std::map<std::string, std::string> argtype;

class Cylinder;
class Table1D;

/**
  A tabular potential for interactions in one cylindrical dimension.
  The table file is formatted as follows:

  As described in ModelTableCart1D,
  the first line should be "site_types=" followed by a comma-separated list of
  the labels of each site type.
  The remaining lines are the individual tables for each of the site types.
  For example, "site_types=O,H" will then contain a line for the table of "O"
  sites, then another line for the table of "H" sites.
  Each table is given in a single line by a number of space-separated values
  that represent the interaction energy spanning the range of the nearest
  distance, r, from the cylindrical surface, scaled by the radius, R, with the
  range z = r/R = [0, 1] and linearly-spaced.
  For example, if the line contains 5 values, then z=r/R=[0, 0.25, 0.5, 0.75, 1]
  and the potential will be interpolated with
  Table1D::forward_difference_interpolation.
  Note that r is not the radial distance.
  Instead, r=0 at the radius (e.g., the inner surface of the cylinder),
  and r=radius at the center.
 */
class ModelTableCylinder1D : public ModelOneBody {
 public:
  //@{
  /** @name Arguments
    - table_file: file name for the table.
    - Cylinder arguments.
   */
  explicit ModelTableCylinder1D(argtype args = argtype());
  explicit ModelTableCylinder1D(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  void precompute(const Configuration& config) override;

  double energy(
    const Position& wrapped_site,
    const Site& site,
    const Configuration& config,
    const ModelParams& model_params) override;

  void serialize(std::ostream& ostr) const override;
  std::shared_ptr<Model> create(std::istream& istr) const override {
    return std::make_shared<ModelTableCylinder1D>(istr); }
  std::shared_ptr<Model> create(argtype * args) const override {
    return std::make_shared<ModelTableCylinder1D>(args); }
  explicit ModelTableCylinder1D(std::istream& istr);
  virtual ~ModelTableCylinder1D();

  //@}
 private:
  std::string table_file_;
  std::unique_ptr<Cylinder> cylinder_;
  std::vector<std::unique_ptr<Table1D> > tables_;
};

}  // namespace feasst

#endif  // FEASST_CONFINEMENT_MODEL_TABLE_CYLINDRICAL_H_
