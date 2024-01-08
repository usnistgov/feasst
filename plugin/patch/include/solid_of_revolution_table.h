
#ifndef FEASST_PATCH_SOLID_OF_REVOLUTION_TABLE_H_
#define FEASST_PATCH_SOLID_OF_REVOLUTION_TABLE_H_

#include <memory>
#include "utils/include/arguments.h"
#include "math/include/table.h"
#include "math/include/matrix.h"
#include "math/include/euler.h"
#include "system/include/visit_model.h"

namespace feasst {

/**
  Model solids of revolution using a tabular potential that is precomputed
  and read from a file.
  An additional assumption is that the solids of revolution are not only
  rotationally symmetric about the axis of revolution, but also symmetric
  about the plane perpendicular to the axis of revolution that passes through
  the center of the particle.
  This class is implemented in a similar fashion as VisitModelInnerAniso,
  except that directors are used for orientations and not euler angles.

  The relative position and orientation of two solids of revolution is given by
  given by three orientational degrees of freedom and one radial distance, r,
  between the centers of the two bodies.

  The three orientational degrees of freedom (i.e., three angles) are defined by
  two unit vectors which represent the orientation of the axis of symmetry of
  each shape (e.g., i2 and j2), and the vector connecting the centers of the two
  shapes, (rij = ri - rj).

       i2            obtain dihedral    i2
      /  thetai          project        |
     i <--- rij ---- j    >>>>>         ij  -psi
           thetaj   /    on rij          \
                  j2                      j2

  This implementation computes -cos(psi), where normally you would expect it
  to be psi.
  The tabular potential uses the -cos(psi) convention(cos(pi-a)=-cos(a)),
  and the nj=j2cross(-uij) leads to -cos(psi).
  Use the right hand rule to see that psi was computed with the negative.
  See: http://math.stackexchange.com/questions/47059/how-do-i-calculate-a-dihedral-angle-given-cartesian-coordinates

  r = |rij|
  uij = rij/r
  cos thetai = i2 .dot. (-uij)
  cos thetaj = j2 .dot. ( uij)
  ni = i2 .cross. uij
  nj = j2 .cross. (-uij)
  cos psi = ni .dot. nj

  Bounds:
  theta(i,j) [0, pi/2] -> symmetry, if theta > pi/2, theta=pi-theta
  psi [0, pi] -> symmetry, if psi > pi, psi = 2pi-psi
  costheta(i,j) ~ [0, 1]
  cospsi ~ [-1, 1]

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

  1. "num_orientations_per_half_pi [value]" is used to determine the resolution of all 5
     orientational degrees of freedom.

  2. "gamma [value]" is the expononent for the definition of z.
     But, if gamma==0, then use a square well potential.

  3. "delta [value]" is the distance between the hard particle inner cutoff and the
     outer cutoff when the potential goes to zero (e.g., delta = r_c - r_h).

  4. "num_z [value]" is the number of energy values along the z parameter.

  5. "smoothing_distance [value]" is the distance to linearly interpolate the
     energy to zero. Specifically, r_c - smoothing_distance is the z=1 value.
     If smoothing_distance is <= 0, then it is not used at all.

  6. The remaining lines are for each unique set of the three angles separated by spaces.
     The outer loop is cos(theta1), with range [0, pi/2].
     The next loop is cos(theta2) in [0, 1], then cos(psi) in [-1, 0].
     The number of lines should be (2k+1)(k+1)^2,
     where k=num_orientations_per_half_pi.
     Each line then reports the contact distance, r_h, followed by num_z values
     of the potential energy uniformly in the z range of [0, 1].

  Coordinate axes and references:

    z
    |
    |__ x    the reference is for the axis of rotation to point along the z-axis
   /         which is how one defines the rotation matrices for the theta_i/j and psi
  y

  In order to obtain the orientation given thetai, thetaj, psi and d:
  -place two particles on the origin with axes of rotation along the z-axis
  -rotate the first particle by an angle phi = pi/2 - thetai about the y axis (toward the x-axis)
  -rotate the second particle by an angle phi = pi/2- thetaj about the y axis (toward the negative x-axis)
  -rotate the second particle by an angle psi about the x-axis
  -translate the second particle by a distance, d along the x-axis

  If using superballs, only the z-exponent, a3 can be different from a1=a2.
  In this case, the reference position has the axis of revolution along the z-axis.
  Therefore, when thetai = PI/2, the particle is in the 'reference' orientation.
 */
class SolidOfRevolutionTable : public VisitModelInner {
 public:
  //@{
  /** @name Arguments
    - table_file: table file with format described above.
    - ignore_energy: do not read the energy table (default: false).
   */
  explicit SolidOfRevolutionTable(argtype args = argtype());
  explicit SolidOfRevolutionTable(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  /// Return true if there is an energy table.
  bool is_energy_table(const std::vector<std::vector<Table4D> >& energy) const;

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
    return std::make_shared<SolidOfRevolutionTable>(istr); }
  std::shared_ptr<VisitModelInner> create(argtype * args) const override {
    return std::make_shared<SolidOfRevolutionTable>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit SolidOfRevolutionTable(std::istream& istr);
  virtual ~SolidOfRevolutionTable() {}

  //@}
 private:
  int director_index_ = -1;
  std::string table_file_;
  bool ignore_energy_;
  std::vector<int> site_types_;
  std::vector<int> t2index_;
  std::vector<std::vector<double> > gamma_;
  std::vector<std::vector<double> > delta_;
  std::vector<std::vector<double> > smoothing_distance_;
  Position dir1_pos_, dir2_pos_;

  void read_table_(const std::string table_file, const bool ignore_energy, Configuration * config);
};

inline std::shared_ptr<SolidOfRevolutionTable> MakeSolidOfRevolutionTable(
    argtype args = argtype()) {
  return std::make_shared<SolidOfRevolutionTable>(args);
}

}  // namespace feasst

#endif  // FEASST_PATCH_SOLID_OF_REVOLUTION_TABLE_H_
