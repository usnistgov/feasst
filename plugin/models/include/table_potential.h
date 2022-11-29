
#ifndef FEASST_MODELS_TABLE_POTENTIAL_H_
#define FEASST_MODELS_TABLE_POTENTIAL_H_

#include <memory>
#include "utils/include/arguments.h"
#include "math/include/table.h"
#include "system/include/visit_model.h"

namespace feasst {

/**
  Represent isotropic interactions using a tabular potential that is precomputed
  and read from a file.

  For an isotropic site, the radial distance, r, between the centers is scaled as
  \f$z=[(r^\gamma - r_h^\gamma)/(r_c^\gamma - r_h^\gamma)]\f$
  which varies from z=0 at the hard contact distance, r_h, to z=1 at the cutoff, r_c.
  A stretching parameter, \f$\gamma\f$, increases the resolution at shorter
  distances when negative.
  This implementation is currently hard coded for \f$\gamma=-2\f$ which avoids sqrt.

  Interpolation uses Table1D::forward_difference_interpolation.

  The format of the tabular potential stored in a file is as follows.

  The first line should be 'site_types' followed by the number of site types
  and then the identity of each of those site types in order of the tables
  given below. (e.g., "site_types n i" where n is the number of site types and
  each following number is the type of each anisotropic site.)

  There is a table for each unique pair of site types.
  For example, if the first line is as follows:
  "site_types 2 1 7"
  then the first table will be for interactions of site type 1 with
  site type 1 (i.e., 1-1), the second table will be for 1-7 interactions,
  and the third table will be for 7-7 interactions.

  For each pair of site types, i <= j, a table is given by the following lines.

  1. "inner [value]" is the inner hard cutoff distance.
  2. "num_z [value]" is the number of energy values along the z parameter.
  3. The last line is num_z space-separated values of the potential
     energy uniformly in the z range of [0, 1].
 */
class TablePotential : public VisitModelInner {
 public:
  /**
    args:
    - table_file: table file with format described above.
   */
  explicit TablePotential(argtype args = argtype());
  explicit TablePotential(argtype * args);

  const std::vector<std::vector<Table1D> >& energy_table() const { return energy_table_; }

  void precompute(Configuration * config) override;
  void compute(
    const int part1_index,
    const int site1_index,
    const int part2_index,
    const int site2_index,
    const Configuration * config,
    const ModelParams& model_params,
    ModelTwoBody * model,
    const bool is_old_config,
    Position * relative,
    Position * pbc) override;

  std::shared_ptr<VisitModelInner> create(std::istream& istr) const override {
    return std::make_shared<TablePotential>(istr); }
  std::shared_ptr<VisitModelInner> create(argtype * args) const override {
    return std::make_shared<TablePotential>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit TablePotential(std::istream& istr);
  virtual ~TablePotential() {}

 private:
  std::vector<std::vector<double> > inner_, inner_g_;
  std::vector<std::vector<double> > cutoff_g_;
  std::vector<int> site_types_;
  //std::vector<std::vector<double> > gamma_;
  std::vector<std::vector<Table1D> > energy_table_;
  const double gamma_ = -2;

  void read_table_(const std::string table_file);
};

inline std::shared_ptr<TablePotential> MakeTablePotential(
    argtype args = argtype()) {
  return std::make_shared<TablePotential>(args);
}

}  // namespace feasst

#endif  // FEASST_MODELS_TABLE_POTENTIAL_H_
