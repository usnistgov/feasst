
#ifndef FEASST_MONTE_CARLO_WRITE_MODEL_PARAMS_H_
#define FEASST_MONTE_CARLO_WRITE_MODEL_PARAMS_H_

#include <vector>
#include <memory>
#include "utils/include/arguments.h"
#include "monte_carlo/include/action.h"

namespace feasst {

/**
  Write ModelParams to file.
 */
class WriteModelParams : public Action {
 public:
  //@{
  /** @name Arguments
    - output_file: name of file to write.
    - potential_index: index of potential.
      If -1, use Configuration (default: -1)
    - reference_index: index of reference potential (default: -1).
   */
  explicit WriteModelParams(argtype args = argtype());
  explicit WriteModelParams(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  void run(MonteCarlo * mc) override;
  std::shared_ptr<Action> create(std::istream& istr) const override {
    return std::make_shared<WriteModelParams>(istr); }
  std::shared_ptr<Action> create(argtype * args) const override {
    return std::make_shared<WriteModelParams>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit WriteModelParams(std::istream& istr);
  virtual ~WriteModelParams() {}

  //@}
 private:
  int potential_index_;
  int reference_index_;
  std::string output_file_;
};

inline std::shared_ptr<WriteModelParams> MakeWriteModelParams(
    argtype args = argtype()) {
  return std::make_shared<WriteModelParams>(args);
}

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_WRITE_MODEL_PARAMS_H_
