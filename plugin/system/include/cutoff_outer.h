
#ifndef FEASST_SYSTEM_CUTOFF_OUTER_H_
#define FEASST_SYSTEM_CUTOFF_OUTER_H_

#include <memory>
#include "configuration/include/model_params.h"

namespace feasst {

class CutoffOuter : public ModelParam {
 public:
  CutoffOuter() : ModelParam() { class_name_ = "cutoff_outer"; }
  std::shared_ptr<ModelParam> create(std::istream& istr) const override {
    return std::make_shared<CutoffOuter>(istr); }
  void serialize(std::ostream& ostr) const override;
  explicit CutoffOuter(std::istream& istr);
  virtual ~CutoffOuter() {}
};

}  // namespace feasst

#endif  // FEASST_SYSTEM_CUTOFF_OUTER_H_
