
#ifndef FEASST_CHAIN_PERTURB_DISTANCE_ANGLE_CONNECTOR_H_
#define FEASST_CHAIN_PERTURB_DISTANCE_ANGLE_CONNECTOR_H_

#include "math/include/euler.h"
#include "monte_carlo/include/perturb_distance_angle.h"

namespace feasst {

/**
  Place an anisotropic (rigid body) site according to the angle and bond,
  as in PerturbDistanceAngle, but also orient the anisotropic site
  according to the position of the first anchor, which is assumed to be
  a connector.
 */
class PerturbDistanceAngleConnector : public PerturbDistanceAngle {
 public:
  explicit PerturbDistanceAngleConnector(argtype args = argtype());
  explicit PerturbDistanceAngleConnector(argtype * args);

  void move_once(const bool is_position_held,
    System * system,
    TrialSelect * select,
    Random * random,
    double * bond_energy) override;

  // serialize
  std::shared_ptr<Perturb> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit PerturbDistanceAngleConnector(std::istream& istr);
  virtual ~PerturbDistanceAngleConnector() {}

 protected:
  void serialize_perturb_distance_angle_connector_(std::ostream& ostr) const;

 private:
  // temporary and not serialized
  RotationMatrix rot_mat_;
  Matrix tmp_mat_;
  Position tmp_pos_;
  Euler euler_;
};

inline std::shared_ptr<PerturbDistanceAngleConnector> MakePerturbDistanceAngleConnector(
    argtype args = argtype()) {
  return std::make_shared<PerturbDistanceAngleConnector>(args);
}

}  // namespace feasst

#endif  // FEASST_CHAIN_PERTURB_DISTANCE_ANGLE_CONNECTOR_H_
