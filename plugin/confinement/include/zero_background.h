
#ifndef FEASST_CONFINEMENT_ZERO_BACKGROUND_H_
#define FEASST_CONFINEMENT_ZERO_BACKGROUND_H_

#include <memory>
#include <vector>
#include "monte_carlo/include/action.h"

namespace feasst {

/**
  Add the negative of the current energy as Background so that the total energy
  is now zero.
 */
class ZeroBackground : public Action {
 public:
  //@{
  /** @name Arguments
    - config: name of Configuration (default: 0).
   */
  explicit ZeroBackground(argtype args = argtype());
  explicit ZeroBackground(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  void run(MonteCarlo * mc) override;
  std::shared_ptr<Action> create(std::istream& istr) const override {
    return std::make_shared<ZeroBackground>(istr); }
  std::shared_ptr<Action> create(argtype * args) const override {
    return std::make_shared<ZeroBackground>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit ZeroBackground(std::istream& istr);
  virtual ~ZeroBackground() {}

  //@}
 private:
  std::string config_;
  int configuration_index_;
};

inline std::shared_ptr<ZeroBackground> MakeZeroBackground(
    argtype args = argtype()) {
  return std::make_shared<ZeroBackground>(args);
}

}  // namespace feasst

#endif  // FEASST_CONFINEMENT_ZERO_BACKGROUND_H_
