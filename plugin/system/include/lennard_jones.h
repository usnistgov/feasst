
#ifndef FEASST_SYSTEM_LENNARD_JONES_H_
#define FEASST_SYSTEM_LENNARD_JONES_H_

#include <map>
#include <string>
#include <memory>
#include "system/include/model_two_body.h"

namespace feasst {

typedef std::map<std::string, std::string> argtype;

/**
  The Lennard-Jones potential is given by

  \f$ U_{LJ} = 4\epsilon [ (\sigma/r)^{2\alpha} - (\sigma/r)^\alpha ] \f$,

  where \f$\alpha=6\f$ is assumed by this class for optimization.

 */
class LennardJones : public ModelTwoBody {
 public:
  //@{
  /** @name Arguments
    - hard_sphere_threshold: when r < threshold*sigma, return NEAR_INFINITY
      (default: 0.2).
   */
  explicit LennardJones(argtype args = argtype());
  explicit LennardJones(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  double energy(
      const double squared_distance,
      const int type1,
      const int type2,
      const ModelParams& model_params) override;

  /// Return the threshold for hard sphere interaction.
  double hard_sphere_threshold() const;
  const double& hard_sphere_threshold_sq() const {
    return hard_sphere_threshold_sq_; }

  std::shared_ptr<Model> create(std::istream& istr) const override {
    return std::make_shared<LennardJones>(istr); }
  std::shared_ptr<Model> create(argtype * args) const override {
    return std::make_shared<LennardJones>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit LennardJones(std::istream& istr);
  virtual ~LennardJones() {}

  //@}
 protected:
  void serialize_lennard_jones_(std::ostream& ostr) const;

 private:
  double hard_sphere_threshold_sq_;
};

inline std::shared_ptr<LennardJones> MakeLennardJones(
    argtype args = argtype()) {
  return std::make_shared<LennardJones>(args);
}

}  // namespace feasst

#endif  // FEASST_SYSTEM_LENNARD_JONES_H_
