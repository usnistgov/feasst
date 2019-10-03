
#ifndef FEASST_MONTE_CARLO_PERTURB_DISTANCE_ANGLE_H_
#define FEASST_MONTE_CARLO_PERTURB_DISTANCE_ANGLE_H_

#include "monte_carlo/include/perturb_distance.h"

namespace feasst {

/// Put first site in selection, i, in a sphere about the first site in anchor,
///  j, and at an angle i,j,k (vertex: j) about the second site in anchor, j.
class PerturbDistanceAngle : public PerturbDistance {
 public:
  PerturbDistanceAngle(const argtype& args = argtype());

  void precompute(TrialSelect * select, System * system) override {
    PerturbDistance::precompute(select, system);
    angle_ = select->property("theta0");
    origin_.set_to_origin(system->configuration().dimension());
    rjk_ = origin_;
  }

  void move(System * system,
      TrialSelect * select,
      Random * random) override {
    SelectList * mobile = select->get_mobile();
    Position * site = mobile->get_site_position(0, 0);
    DEBUG("mobile " << mobile->str());
    DEBUG("old pos " << site->str());

    // set site to the vector |r_j - r_k| and store this unit vector
    const Position& rj = select->anchor_position(0, 0, system);
    const Position& rk = select->anchor_position(0, 1, system);
    *site = rj;
    site->subtract(rk);
    rjk_ = *site;
    DEBUG("rjk " << rjk_.str());
    //rjk_.normalize();

    // rotate site by (PI-theta0) about vector orthogonal to r_jk
    orthogonal_jk_.orthogonal(*site);
    DEBUG("ortho " << orthogonal_jk_.str());
    rot_mat_.axis_angle(orthogonal_jk_, PI - angle_);
    DEBUG("site == rj: " << site->str());
    rot_mat_.rotate(origin_, site);
    DEBUG("site rotated to angle: " << site->str());

    // randomly spin site about rjk.
    rot_mat_.axis_angle(rjk_, 2*PI*random->uniform());
    rot_mat_.rotate(origin_, site);

    site->add(rj);  // return frame of reference

    DEBUG("new pos " << site->str());
    system->get_configuration()->update_positions(select->mobile());
  }

  std::shared_ptr<Perturb> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit PerturbDistanceAngle(std::istream& istr);
  virtual ~PerturbDistanceAngle() {}

 private:
  double angle_ = 0.;

  // temporary
  Position rjk_;
  Position orthogonal_jk_;
  Position origin_;
  RotationMatrix rot_mat_;
};

inline std::shared_ptr<PerturbDistanceAngle> MakePerturbDistanceAngle(const argtype& args = argtype()) {
  return std::make_shared<PerturbDistanceAngle>(args);
}

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_PERTURB_DISTANCE_ANGLE_H_
