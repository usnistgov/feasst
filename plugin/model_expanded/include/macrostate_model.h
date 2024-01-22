
#ifndef FEASST_MODEL_EXPANDED_MACROSTATE_MODEL_H_
#define FEASST_MODEL_EXPANDED_MACROSTATE_MODEL_H_

#include "flat_histogram/include/macrostate.h"
#include "model_expanded/include/constrain_model_index.h"

namespace feasst {

/**
  Defines the macrostate as ModelExpanded::model_index.
 */
class MacrostateModel : public Macrostate {
 public:
  //@{
  /** @name Arguments
    - ConstrainModelIndex arguments.
    - Macrostate arguments.
   */
  explicit MacrostateModel(argtype args = argtype());
  explicit MacrostateModel(argtype * args) :
    MacrostateModel(Histogram(args), args) {}
  //@}
  /** @name Public Functions
   */
  //@{
  MacrostateModel(const Histogram& histogram, argtype args = argtype());
  MacrostateModel(const Histogram& histogram, argtype * args);
  double value(const System& system,
    const Criteria& criteria,
    const Acceptance& acceptance) const override;
  std::shared_ptr<Macrostate> create(std::istream& istr) const override;
  std::shared_ptr<Macrostate> create(argtype * args) const override;
  void serialize(std::ostream& ostr) const override;
  MacrostateModel(std::istream& istr);
  virtual ~MacrostateModel() {}
  //@}

 private:
  ConstrainModelIndex constraint_;
};

inline std::shared_ptr<MacrostateModel> MakeMacrostateModel(
    const Histogram& histogram, argtype args = argtype()) {
  return std::make_shared<MacrostateModel>(histogram, args);
}

inline std::shared_ptr<MacrostateModel> MakeMacrostateModel(
    argtype args = argtype()) {
  return std::make_shared<MacrostateModel>(args);
}

}  // namespace feasst

#endif  // FEASST_MODEL_EXPANDED_MACROSTATE_MODEL_H_
