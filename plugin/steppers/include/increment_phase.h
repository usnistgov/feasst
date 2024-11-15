
#ifndef FEASST_STEPPERS_INCREMENT_PHASE_H_
#define FEASST_STEPPERS_INCREMENT_PHASE_H_

#include "monte_carlo/include/modify.h"

namespace feasst {

/**
  Increment the simulation phase in Criteria when input conditions are reached.
 */
class IncrementPhase : public ModifyUpdateOnly {
 public:
  //@{
  /** @name Arguments
    - num_trials: increment the simulation phase when this many trials have
      been attempted. Do nothing if -1 (default: -1).
  */
  explicit IncrementPhase(argtype args = argtype());
  explicit IncrementPhase(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  void update(MonteCarlo * mc) override;

  std::string class_name() const override { return std::string("IncrementPhase"); }

  // serialize
  void serialize(std::ostream& ostr) const override;
  std::shared_ptr<Modify> create(std::istream& istr) const override {
    return std::make_shared<IncrementPhase>(istr); }
  std::shared_ptr<Modify> create(argtype * args) const override {
    return std::make_shared<IncrementPhase>(args); }
  explicit IncrementPhase(std::istream& istr);

 private:
  int num_trials_;
};

inline std::shared_ptr<IncrementPhase> MakeIncrementPhase(
    argtype args = argtype()) {
  return std::make_shared<IncrementPhase>(args);
}

}  // namespace feasst

#endif  // FEASST_STEPPERS_INCREMENT_PHASE_H_
