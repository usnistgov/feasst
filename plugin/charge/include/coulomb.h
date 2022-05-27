
#ifndef FEASST_CHARGE_COULOMB_H_
#define FEASST_CHARGE_COULOMB_H_

#include "system/include/model_two_body.h"

namespace feasst {

/**
  Compute energy between two point charges, \f$q_i\f$ and \f$q_j\f$ with a
  plain Coulomb potential.

  \f$U = q_i q_j \chi/r\f$

  where \f$r\f$ is the separation distance,
  and \f$\chi\f$ is the charge conversion factor assuming the following units:
    1. length: Angstroms
    2. energy: kJ/mol
    3. charge: elementary

  Avoid Coulomb explosion by returning a large number when \f$r\f$ is near zero.
 */
class Coulomb : public ModelTwoBody {
 public:
  Coulomb() { class_name_ = "Coulomb"; }

  double energy(
      const double squared_distance,
      const int type1,
      const int type2,
      const ModelParams& model_params) override;

  void precompute(const ModelParams& existing) override;

  std::shared_ptr<Model> create(std::istream& istr) const override {
    return std::make_shared<Coulomb>(istr); }
  std::shared_ptr<Model> create(argtype * args) const override {
    return std::make_shared<Coulomb>(); }
  void serialize(std::ostream& ostr) const override;
  explicit Coulomb(std::istream& istr);
  virtual ~Coulomb() {}

 private:
  double conversion_factor_;
};

inline std::shared_ptr<Coulomb> MakeCoulomb() {
  return std::make_shared<Coulomb>();
}

}  // namespace feasst

#endif  // FEASST_CHARGE_COULOMB_H_
