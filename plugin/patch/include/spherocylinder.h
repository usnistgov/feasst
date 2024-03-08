
#ifndef FEASST_PATCH_SPHEROCYLINDER_H_
#define FEASST_PATCH_SPHEROCYLINDER_H_

#include "utils/include/arguments.h"
#include "system/include/visit_model.h"
#include "patch/include/patch_angle.h"

namespace feasst {

/**
  Implementation to find the nearest distance between two spherocylinders
  as described in

  http://dx.doi.org/10.1016/0097-8485(94)80023-5

  and as implemented in

  http://www.sklogwiki.org/SklogWiki/index.php/Source_code_for_the_minimum_distance_between_two_rods#cite_note-1

  but more specifically with the C variation that uses symmetry:

  http://www.sklogwiki.org/SklogWiki/index.php/Rev._source_code_for_the_minimum_distance_between_two_rods_in_C

  A spherocylinder is defined by at two sites that are bonded together.
  The first is the center of the spherocylinder.
  The center should have a cutoff that is the maximum possible interaction distance.
  Thus, the cutoff for the center should be the maximum length + the maximum cutoff of the remaining sites.
  The remaining sites are the orientations of the spherocylinder (e.g., director).
  The cutoff should be the cutoff applied to the shortest distance between spherocylinders
  The sigma and epsilons also play a role (e.g., for HardSphere, SquareWell).
  The second site also contains the length of the spherocylinder.
 */
class Spherocylinder : public VisitModelInner {
 public:
  explicit Spherocylinder(argtype args = argtype());
  explicit Spherocylinder(argtype * args);
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
    Position * pbc,
    const double weight = 1.) override;

  const SpherocylinderLength& length() const { return length_; }

  std::shared_ptr<VisitModelInner> create(std::istream& istr) const override {
    return std::make_shared<Spherocylinder>(istr); }
  std::shared_ptr<VisitModelInner> create(argtype * args) const override {
    return std::make_shared<Spherocylinder>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit Spherocylinder(std::istream& istr);
  virtual ~Spherocylinder() {}

 private:
  SpherocylinderLength length_;
  int director_index_ = -1;
  Position dir1_pos_, dir2_pos_;

  // temporary and not serialized
  double xla_, xmu_;
  double calc_sph_sq_dist_vega_(const double rr, const double rw1,
    const double rw2, const double w1w2, const double lh1, const double lh2);
};

inline std::shared_ptr<Spherocylinder> MakeSpherocylinder(
    argtype args = argtype()) {
  return std::make_shared<Spherocylinder>(args);
}

}  // namespace feasst

#endif  // FEASST_PATCH_SPHEROCYLINDER_H_
