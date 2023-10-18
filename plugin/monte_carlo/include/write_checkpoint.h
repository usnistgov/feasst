
#ifndef FEASST_MONTE_CARLO_WRITE_CHECKPOINT_H_
#define FEASST_MONTE_CARLO_WRITE_CHECKPOINT_H_

#include <vector>
#include <memory>
#include "utils/include/arguments.h"
#include "monte_carlo/include/action.h"

namespace feasst {

/**
  Write a Checkpoint.
 */
class WriteCheckpoint : public Action {
 public:
  explicit WriteCheckpoint(argtype args = argtype());
  explicit WriteCheckpoint(argtype * args);
  void run(MonteCarlo * mc) override;
  std::shared_ptr<Action> create(std::istream& istr) const override {
    return std::make_shared<WriteCheckpoint>(istr); }
  std::shared_ptr<Action> create(argtype * args) const override {
    return std::make_shared<WriteCheckpoint>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit WriteCheckpoint(std::istream& istr);
  virtual ~WriteCheckpoint() {}
};

inline std::shared_ptr<WriteCheckpoint> MakeWriteCheckpoint(argtype args = argtype()) {
  return std::make_shared<WriteCheckpoint>(args);
}

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_WRITE_CHECKPOINT_H_
