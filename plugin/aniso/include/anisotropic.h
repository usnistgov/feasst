
#ifndef FEASST_ANISO_ANISOTROPIC_H_
#define FEASST_ANISO_ANISOTROPIC_H_

#include "configuration/include/model_params.h"

namespace feasst {

/**
  The anisotropic parameter in LMP-like data file Site Properties signals that
  a site has orientation and must track its Euler angles.
 */
class Anisotropic : public ModelParam {
 public:
  Anisotropic() : ModelParam() { class_name_ = "anisotropic"; }
  std::shared_ptr<ModelParam> create(std::istream& istr) const override {
    return std::make_shared<Anisotropic>(istr); }
  void serialize(std::ostream& ostr) const override;
  explicit Anisotropic(std::istream& istr);
  virtual ~Anisotropic() {}
};

}  // namespace feasst

#endif  // FEASST_ANISO_ANISOTROPIC_H_
