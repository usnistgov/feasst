
#ifndef FEASST_MODELS_TABLE_POTENTIAL_H_
#define FEASST_MODELS_TABLE_POTENTIAL_H_

#include <memory>
#include "math/include/table.h"
#include "system/include/model_two_body.h"

namespace feasst {

typedef std::map<std::string, std::string> argtype;

/**
  Represent isotropic interactions using a tabular potential that is precomputed
  and read from a file.

  For an isotropic site, the radial distance, r, between the centers is scaled as
  \f$z=[(r^\gamma - r_h^\gamma)/(r_c^\gamma - r_h^\gamma)]\f$
  which varies from z=0 at the hard contact distance, \f$r_h\f$, to \f$z=1\f$ at
  the cutoff, \f$r_c\f$.
  A stretching parameter, \f$\gamma\f$, increases the resolution at shorter
  distances when negative.

  The default value of \f$\gamma=-2\f$ avoids the expensive sqrt operation, but
  can also be problematic if \f$r_h \approx 0\f$.
  For \f$r_h \approx 0\f$, \f$\gamma=1\f$ may be more appropriate.

  Interpolation uses Table1D::forward_difference_interpolation.

  The format of the tabular potential stored in a file is as follows.

  The first line should be 'site_types' followed by the number of site types
  and then the name of each of those site types in order of the tables
  given below. (e.g., "site_types 1 name" where n is the number of site types and
  each following name is the site type name of each site.)

  There is a table for each unique pair of site types.
  For example, if the first line is as follows:
  "site_types 2 O H"
  then the first table will be for interactions of site type O with
  site type O (i.e., O-O), the second table will be for O-H interactions,
  and the third table will be for H-H interactions.

  For each pair of site types, a table is given by the following lines.

  1. "gamma [value]" is the optional stretching exponential.
     If this line is not provided, then the default value of -2 is used.
  2. "inner [value]" is the inner hard cutoff distance.
  3. "num_z [value]" is the number of energy values along the z parameter.
  4. The last line is num_z space-separated values of the potential
     energy uniformly in the z range of [0, 1].
 */
class TablePotential : public ModelTwoBody {
 public:
  //@{
  /** @name Arguments
    - table_file: table file with format described above.
   */
  explicit TablePotential(argtype args = argtype());
  explicit TablePotential(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  const std::vector<std::vector<Table1D> >& energy_table() const { return energy_table_; }

  void precompute(const Configuration& config) override;
  double energy(
    const double squared_distance,
    const int type1,
    const int type2,
    const ModelParams& model_params) override;

  std::shared_ptr<Model> create(std::istream& istr) const override {
    return std::make_shared<TablePotential>(istr); }
  std::shared_ptr<Model> create(argtype * args) const override {
    return std::make_shared<TablePotential>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit TablePotential(std::istream& istr);
  virtual ~TablePotential() {}

  //@}
 private:
  std::vector<std::vector<double> > inner_, inner_g_;
  std::vector<std::vector<double> > cutoff_g_;
  std::vector<int> site_types_;
  std::vector<std::string> site_type_names_;
  std::vector<int> t2index_;
  std::vector<std::vector<double> > gamma_;
  std::vector<std::vector<Table1D> > energy_table_;

  void read_table_(const std::string table_file);
  void read_inner_(const int itype, const int jtype, const std::string& description, const double value, std::ifstream& file);
};

inline std::shared_ptr<TablePotential> MakeTablePotential(
    argtype args = argtype()) {
  return std::make_shared<TablePotential>(args);
}

}  // namespace feasst

#endif  // FEASST_MODELS_TABLE_POTENTIAL_H_
