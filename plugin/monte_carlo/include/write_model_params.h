
#ifndef FEASST_MONTE_CARLO_WRITE_MODEL_PARAMS_H_
#define FEASST_MONTE_CARLO_WRITE_MODEL_PARAMS_H_

#include <map>
#include <memory>
#include <string>
#include <vector>
#include "monte_carlo/include/action.h"

namespace feasst {

typedef std::map<std::string, std::string> argtype;

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
    - config: name of Configuration (default: 0).
    - ref: name of RefPotential. If empty, use Potential (default: empty).
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
  std::string config_;
  std::string ref_;
  std::string output_file_;
};

inline std::shared_ptr<WriteModelParams> MakeWriteModelParams(
    argtype args = argtype()) {
  return std::make_shared<WriteModelParams>(args);
}

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_WRITE_MODEL_PARAMS_H_
