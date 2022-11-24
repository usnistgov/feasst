
#ifndef FEASST_MODELS_TWO_BODY_TABLE_H_
#define FEASST_MODELS_TWO_BODY_TABLE_H_

#include <string>
#include <memory>
#include "system/include/model_two_body.h"

namespace feasst {

/**
  Use this as the Model for TablePotential.
 */
class TwoBodyTable : public ModelTwoBody {
 public:
  TwoBodyTable() { class_name_ = "TwoBodyTable"; }

  double energy(
      const double squared_distance,
      const int type1,
      const int type2,
      const ModelParams& model_params) override;

  std::shared_ptr<Model> create(std::istream& istr) const override {
    return std::make_shared<TwoBodyTable>(istr); }
  std::shared_ptr<Model> create(argtype * args) const override {
    return std::make_shared<TwoBodyTable>(); }
  void serialize(std::ostream& ostr) const override;
  explicit TwoBodyTable(std::istream& istr);
  virtual ~TwoBodyTable() {}

 protected:
  void serialize_ideal_gas_(std::ostream& ostr) const;
};

inline std::shared_ptr<TwoBodyTable> MakeTwoBodyTable() {
  return std::make_shared<TwoBodyTable>();
}

}  // namespace feasst

#endif  // FEASST_MODELS_TWO_BODY_TABLE_H_
