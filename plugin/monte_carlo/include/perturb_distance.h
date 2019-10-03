
#ifndef FEASST_MONTE_CARLO_PERTURB_DISTANCE_H_
#define FEASST_MONTE_CARLO_PERTURB_DISTANCE_H_

#include "monte_carlo/include/perturb_move.h"

namespace feasst {

/// Put first site in selection in a sphere about the first site in anchor.
/// HWH: could enable tuning and position w.r.t. previous bond placement
///      for higher acceptance probability
class PerturbDistance : public PerturbMove {
 public:
  PerturbDistance(const argtype& args = argtype());

  void precompute(TrialSelect * select, System * system) override {
    // determine the bond length
    // or input the bond length
    if (select->has_property("bond_length")) {
      distance_ = select->property("bond_length");
    } else {
      WARN("using default distance (typically for reptation): " << distance_);
    }
  }

  double distance() const { return distance_; }

  void move(System * system,
      TrialSelect * select,
      Random * random) override {
    SelectList * mobile = select->get_mobile();
    Position * site = mobile->get_site_position(0, 0);
    DEBUG("mobile " << mobile->str());
    DEBUG("old pos " << site->str());
    random->unit_sphere_surface(site);
    site->multiply(distance_);
    site->add(select->anchor_position(0, 0, system));
    DEBUG("new pos " << site->str());
    system->get_configuration()->update_positions(select->mobile());
  }

  std::shared_ptr<Perturb> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit PerturbDistance(std::istream& istr);
  virtual ~PerturbDistance() {}

 protected:
  void serialize_perturb_distance_(std::ostream& ostr) const;

 private:
  double distance_ = 1.;
};

inline std::shared_ptr<PerturbDistance> MakePerturbDistance(const argtype& args = argtype()) {
  return std::make_shared<PerturbDistance>(args);
}

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_PERTURB_DISTANCE_H_
