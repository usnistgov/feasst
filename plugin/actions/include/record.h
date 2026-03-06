
#ifndef FEASST_MONTE_CARLO_RECORD_H_
#define FEASST_MONTE_CARLO_RECORD_H_

#include <memory>
#include <string>
#include <vector>
#include <cstdint>
#include "monte_carlo/include/action.h"

namespace feasst {

/**
  Make or save a Record and then load it back later.
 */
class Record : public Action {
 public:
  //@{
  /** @name Arguments
    - config: name of Configuration (default: 0)
    - save_positions: The file name to write the positions of the given
      Configuration, if not empty (default: empty).
    - load_positions: The file name to read and load into the given
      Configuration, if not empty (default: empty).
   */
  explicit Record(argtype args = argtype());
  explicit Record(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  void run(MonteCarlo * mc) override;
  std::shared_ptr<Action> create(std::istream& istr) const override {
    return std::make_shared<Record>(istr); }
  std::shared_ptr<Action> create(argtype * args) const override {
    return std::make_shared<Record>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit Record(std::istream& istr);
  virtual ~Record() {}

  //@}
 private:
  std::string config_, save_positions_, load_positions_;
};

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_RECORD_H_
