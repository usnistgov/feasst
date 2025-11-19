
#ifndef FEASST_ANISO_TRAINING_H_
#define FEASST_ANISO_TRAINING_H_

#include <map>
#include <memory>
#include <string>
#include <vector>
#include "monte_carlo/include/action.h"

namespace feasst {

typedef std::map<std::string, std::string> argtype;

/**
 */
class Training : public Action {
 public:
  //@{
  /** @name Arguments
    - training_file: file name of comma-separated values of 3D rigid body degrees of freedom (spherical, euler and energy).
   */
  explicit Training(argtype args = argtype());
  explicit Training(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  void run(MonteCarlo * mc) override;
  std::shared_ptr<Action> create(std::istream& istr) const override {
    return std::make_shared<Training>(istr); }
  std::shared_ptr<Action> create(argtype * args) const override {
    return std::make_shared<Training>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit Training(std::istream& istr);
  virtual ~Training() {}
  //@}

 private:
  std::string training_file_;
};

}  // namespace feasst

#endif  // FEASST_ANISO_TRAINING_H_
