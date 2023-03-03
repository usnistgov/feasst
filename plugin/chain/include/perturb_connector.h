
#ifndef FEASST_CHAIN_PERTURB_CONNECTOR_H_
#define FEASST_CHAIN_PERTURB_CONNECTOR_H_

#include "math/include/matrix.h"
#include "monte_carlo/include/perturb_move.h"

namespace feasst {

/**
  Place a connector site according to the orientation of an anisotropic site.
  Connector sites typically have no iteractions and are used to model the
  flexibility between anisotropic sites and other sites (both isotropic and
  anisotropic).
 */
class PerturbConnector : public PerturbMove {
 public:
  explicit PerturbConnector(argtype args = argtype());
  explicit PerturbConnector(argtype * args);

  void move(const bool is_position_held, System * system, TrialSelect * select,
            Random * random) override;

  // serialize
  std::shared_ptr<Perturb> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit PerturbConnector(std::istream& istr);
  virtual ~PerturbConnector() {}

 protected:
  // temporary and not serialized
  RotationMatrix rot_mat_;  

  void serialize_perturb_connector_(std::ostream& ostr) const;
};

inline std::shared_ptr<PerturbConnector> MakePerturbConnector(
    argtype args = argtype()) {
  return std::make_shared<PerturbConnector>(args);
}

}  // namespace feasst

#endif  // FEASST_CHAIN_PERTURB_CONNECTOR_H_
