
#ifndef FEASST_CONFINEMENT_BACKGROUND_H_
#define FEASST_CONFINEMENT_BACKGROUND_H_

#include <memory>
#include "system/include/visit_model.h"

namespace feasst {

typedef std::map<std::string, std::string> argtype;

/**
  Add a constant background potential that doesn't depend on the configuration.
 */
class Background : public VisitModel {
 public:
  //@{
  /** @name Arguments
    - constant: return this energy.
   */
  explicit Background(argtype args = argtype());
  //@}
  /** @name Public Functions
   */
  //@{
  void compute(
      ModelOneBody * model,
      const ModelParams& model_params,
      const Select& selection,
      Configuration * config,
      const int group_index) override;
  void compute(
      ModelOneBody * model,
      const ModelParams& model_params,
      Configuration * config,
      const int group_index) override;
  std::shared_ptr<VisitModel> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit Background(std::istream& istr);
  ~Background() {}

  //@}
 private:
  double constant_;
};

inline std::shared_ptr<Background> MakeBackground(
    argtype args = argtype()) {
  return std::make_shared<Background>(args);
}

}  // namespace feasst

#endif  // FEASST_CONFINEMENT_BACKGROUND_H_
