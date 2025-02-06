
#ifndef FEASST_MONTE_CARLO_RUN_H_
#define FEASST_MONTE_CARLO_RUN_H_

#include <memory>
#include <string>
#include <vector>
#include "monte_carlo/include/action.h"

namespace feasst {

/**
  Write a file with the name of prefix + sim + suffix.
  Then, write an empty output_file if all sim of a given num_sim files
  are present.
  This can be used to terminate a simulation after all simulations in a group
  are finished.
 */
class WriteFileAndCheck : public Action {
 public:
  //@{
  /** @name Arguments
    - sim: the index of the simulation starting with 0 to num_sims - 1.
    - file_prefix: the beginning of the file name.
    - file_suffix: the end of the file name.
    - sim_start: the index of the first simulation of joint termination group.
    - sim_end: the index of the last simulation of joint termination group.
    - output_file: the name of the file written if all sim files are present.
   */
  explicit WriteFileAndCheck(argtype args = argtype());
  explicit WriteFileAndCheck(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  void run(MonteCarlo * mc) override;
  std::shared_ptr<Action> create(std::istream& istr) const override {
    return std::make_shared<WriteFileAndCheck>(istr); }
  std::shared_ptr<Action> create(argtype * args) const override {
    return std::make_shared<WriteFileAndCheck>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit WriteFileAndCheck(std::istream& istr);
  virtual ~WriteFileAndCheck() {}

  //@}
 private:
   int sim_, sim_start_, sim_end_;
   std::string file_prefix_, file_suffix_, output_file_;
};

inline std::shared_ptr<WriteFileAndCheck> MakeWriteFileAndCheck(argtype args = argtype()) {
  return std::make_shared<WriteFileAndCheck>(args);
}

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_RUN_H_
