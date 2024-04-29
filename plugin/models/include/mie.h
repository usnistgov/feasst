
#ifndef FEASST_MODELS_MIE_H_
#define FEASST_MODELS_MIE_H_

#include "utils/include/arguments.h"
#include "configuration/include/model_params.h"
#include "system/include/model_two_body.h"

namespace feasst {

/**
  The Mie potential, \f$U\f$ is described in
  http://www.sklogwiki.org/SklogWiki/index.php/Mie_potential.

  \f$U=\epsilon\left(\frac{n}{n-m}\right)\left(\frac{n}{m}\right)^{m/(n-m)}\left[\left(\frac{\sigma}{r}\right)^n-\left(\frac{\sigma}{r}\right)^m\right]\f$
 
  Set n and m as the parameters mie_lambda_r and mie_lambda_a, respectively,
  in "Site Properties" particle files, Configuration or Potential arguments.
 */
class Mie : public ModelTwoBody {
 public:
  //@{
  explicit Mie(argtype args = argtype());
  explicit Mie(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  void precompute(const ModelParams& existing) override;
  
  double energy(
      const double squared_distance,
      const int type1,
      const int type2,
      const ModelParams& model_params) override;

  std::shared_ptr<Model> create(std::istream& istr) const override {
    return std::make_shared<Mie>(istr); }
  std::shared_ptr<Model> create(argtype * args) const override {
    return std::make_shared<Mie>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit Mie(std::istream& istr);
  virtual ~Mie() {}

  //@}
 protected:
  void serialize_mie_(std::ostream& ostr) const;

 private:
  int mie_lambda_r_index_ = -1;
  int mie_lambda_a_index_ = -1;
};

inline std::shared_ptr<Mie> MakeMie(argtype args = argtype()) {
  return std::make_shared<Mie>(args);
}

}  // namespace feasst

#endif  // FEASST_MODELS_MIE_H_
