
#ifndef FEASST_SYSTEM_MODEL_EMPTY_H_
#define FEASST_SYSTEM_MODEL_EMPTY_H_

#include <memory>
#include <string>
#include "system/include/model_one_body.h"

namespace feasst {

/// This model throws an exception if used.
class ModelEmpty : public ModelOneBody {
 public:
  ModelEmpty() { class_name_ = "ModelEmpty"; }

  double energy(
      const Position& wrapped_site,
      const Site& site,
      const Configuration& config,
      const ModelParams& model_params) override;

  std::shared_ptr<Model> create(std::istream& istr) const override {
    return std::make_shared<ModelEmpty>(istr);
  }
  void serialize(std::ostream& ostr) const override;
  explicit ModelEmpty(std::istream& istr);
  virtual ~ModelEmpty() {}
};

inline std::shared_ptr<ModelEmpty> MakeModelEmpty() {
  return std::make_shared<ModelEmpty>();
}

}  // namespace feasst

#endif  // FEASST_SYSTEM_MODEL_EMPTY_H_
