
#ifndef FEASST_ANISO_VISIT_MODEL_INNER_TABLE_H_
#define FEASST_ANISO_VISIT_MODEL_INNER_TABLE_H_

#include <memory>
#include "utils/include/arguments.h"
#include "math/include/table.h"
#include "math/include/matrix.h"
#include "math/include/euler.h"
#include "system/include/visit_model.h"

namespace feasst {

/**
  Represent anisotropic sites using a tabular potential that is precomputed
  and read from a file.

  The relative position and orientation of an arbitrary rigid body is given
  by 5 orientational degrees of freedom and the radial distance, r between
  the centers of the two bodies.

  The spherical coordinate azimuthal angle [0, 2pi] and polar angle [0, pi]
  are denoted s1 and s2, and are as described in the Position class.
  Note that there is increased resolution at the poles.
  While something resembling a Fibonacci lattice could more uniformly represent
  a spherical surface, interpolation would be more complicated.

  The Euler angles, [e1, e2, e3], in order, are as described in the Euler class.

  The radial distance, r, is transformed into z which varies from z=0 at the
  hard contact distance, r_h, to z=1 at the cutoff, r_c.
  A stretching parameter, \f$\gamma\f$, increases the resolution at shorter
  distances when negative.

  \f$z=[(r^\gamma - r_h^\gamma)/(r_c^\gamma - r_h^\gamma)]\f$.

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

  1. "num_orientations_per_pi [value]" is used to determine the resolution of all 5
     orientational degrees of freedom.

  2. "gamma [value]" is the expononent for the definition of z.
     But, if gamma==0, then use a square well potential.

  3. "delta [value]" is the distance between the hard particle inner cutoff and the
     outer cutoff when the potential goes to zero (e.g., delta = r_c - r_h).

  4. "num_z [value]" is the number of energy values along the z parameter.

  5. "smoothing_distance [value]" is the distance to linearly interpolate the
     energy to zero. Specifically, r_c - smoothing_distance is the z=1 value.
     If smoothing_distance is <= 0, then it is not used at all.

  6. The remaining lines are for each unique set of the five angles separated by spaces.
     The outer loop is s1, with range [-pi, pi]. If i == j, the range is [0, pi]
     for i-j swap symmetry (e.g., swap i and j if s1 < 0).
     The next loop is s2 in [0, pi], then e1 in [-pi, pi], e2 and [0, pi],
     and finally e3 in [-pi, pi].
     The number of lines should be (k+1)^3(2k+1)^2 for i==j,
     and (k+1)^2(2k+1)^3 lines for i!=j, where k=num_orientations_per_pi.
     Each line then reports the contact distance, r_h, followed by num_z values
     of the potential energy uniformly in the z range of [0, 1].
     If the orientation is duplicate, then the first number is -1 and the second
     is the integer number of the orientation that it is a duplicate of.
 */
class VisitModelInnerTable : public VisitModelInner {
 public:
  //@{
  /** @name Arguments
   */

  /**
    args:
    - table_file: table file with format described above.
    - ignore_energy: do not read the energy table (default: false).
   */
  explicit VisitModelInnerTable(argtype args = argtype());
  explicit VisitModelInnerTable(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  /// Return true if there is an energy table.
  bool is_energy_table(const std::vector<std::vector<Table6D> >& energy) const;

  /**
    Compute the second virial coefficient between anisotropic sites by table
    integration.
    See https://doi.org/10.1063/1.5040252
    For the hard particle interaction, see Eq D5.

    For the Jacobian in integration, see
    https://doi.org/10.1021%2Facs.jpcb.9b06808
    https://doi.org/10.1063/1.432162 (missing alpha integration?)
    https://doi.org/10.1098/rspa.1948.0009
    https://doi.org/10.1063/1.4961541

    args:
    - site_type1: type of first site (default: 0).
    - site_type2: type of second site (default: 0).
    - expand_t: increase orientation points by this factor (default: 1).
    - expand_z: increase translation points by this factor (default: 1).
    - beta: see ThermoParams (default: 1).
   */
  double second_virial_coefficient(const Configuration& config, argtype args = argtype()) const;

// HWH this doesn't work but may be salvaged later.
//  /**
//    Write an xyz file that represents the shape based on the inner cutoff, rh.
//    For a given pair of types, and for each angle, a site is created at rh/2.
//
//    args:
//    - xyz_file: name of file to print xyz representation
//    - type1: pairwise site type (default: 0).
//    - type2: pairwise site type (default: 0).
//    - expand_t: increase orientation points by this factor (default: 1).
//   */
//  void write_surface(argtype args) const;

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
    return std::make_shared<VisitModelInnerTable>(istr); }
  std::shared_ptr<VisitModelInner> create(argtype * args) const override {
    return std::make_shared<VisitModelInnerTable>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit VisitModelInnerTable(std::istream& istr);
  virtual ~VisitModelInnerTable() {}

  //@}
 private:
  int aniso_index_ = -1;
  std::string table_file_;
  bool ignore_energy_;
  std::vector<int> site_types_;
  std::vector<int> t2index_;
  std::vector<std::vector<double> > gamma_;
  std::vector<std::vector<double> > delta_;
  std::vector<std::vector<double> > smoothing_distance_;

  // no serialized optimization variables
  Position pos1_, pos2_, sph_;
  RotationMatrix rot1_, rot2_, rot3_;
  Euler euler_;

  void read_table_(const std::string table_file, const bool ignore_energy, Configuration * config);
};

inline std::shared_ptr<VisitModelInnerTable> MakeVisitModelInnerTable(
    argtype args = argtype()) {
  return std::make_shared<VisitModelInnerTable>(args);
}

}  // namespace feasst

#endif  // FEASST_ANISO_VISIT_MODEL_INNER_TABLE_H_
