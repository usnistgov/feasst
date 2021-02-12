
#ifndef FEASST_MODELS_MIE_H_
#define FEASST_MODELS_MIE_H_

#include "utils/include/arguments.h"
#include "configuration/include/model_params.h"
#include "system/include/model_two_body.h"

namespace feasst {

/**
  The Mie potential, \f$U\f$ is described in
  http://www.sklogwiki.org/SklogWiki/index.php/Mie_potential.

  /f$U=\epsilon\left(\frac{n}{n-m}\right)\left(\frac{n}{m}\right)^{m/(n-m)}\left[\left(\frac{\sigma}{r}\right)^n-\left(\frac{\sigma}{r}\right)^m\right]
 */
class Mie : public ModelTwoBody {
 public:
  /**
    args:
    - n: set the value of \f$n\f$ (default: 12).
    - m: set the value of \f$m\f$ (default: 6).
   */
  Mie(const argtype& args = argtype());

  /// Return the value of n.
  const double& n() const { return n_; }

  /// Return the value of m.
  const double& m() const { return m_; }

  double energy(
      const double squared_distance,
      const int type1,
      const int type2,
      const ModelParams& model_params) override;

  std::shared_ptr<Model> create(std::istream& istr) const override {
    return std::make_shared<Mie>(istr);
  }

  void serialize(std::ostream& ostr) const override;
  explicit Mie(std::istream& istr);
  virtual ~Mie() {}

 protected:
  void serialize_mie_(std::ostream& ostr) const;

 private:
  double n_;
  double m_;
  double prefactor_;
};

inline std::shared_ptr<Mie> MakeMie(
    const argtype& args = argtype()) {
  return std::make_shared<Mie>(args);
}

}  // namespace feasst

#endif  // FEASST_MODELS_MIE_H_
